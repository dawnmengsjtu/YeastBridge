"""Regression tests for the current manuscript release and its portable installer."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import paper_release as paper


class PaperReleaseTests(unittest.TestCase):
    def test_compatibility_entry_has_no_removed_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            command = [sys.executable, str(ROOT / 'scripts/run_esm2joint_chain.py')]
            help_result = subprocess.run(command + ['--help'], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(help_result.returncode, 0, help_result.stderr)
            invalid = subprocess.run(command + ['--threads', '0'], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(invalid.returncode, 2)
            self.assertIn('--threads must be positive', invalid.stderr)

    def test_export_from_unrelated_directory_is_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for script, extra in [('predict.py', []), ('main.py', ['--mode', 'export'])]:
                out = tmp / script / 'results.csv'
                subprocess.run([sys.executable, str(ROOT / script), *extra, '--output', str(out)],
                               cwd=tmp, check=True, capture_output=True, text=True)
                self.assertEqual(out.read_bytes(), (ROOT / 'data/paper/results.csv').read_bytes())
                meta = json.loads(out.with_suffix('.metadata.json').read_text())
                self.assertEqual((meta['rows'], meta['targets'], meta['compounds']), (25, 7, 18))
                self.assertEqual(meta['mode'], 'published-export')
                self.assertEqual(meta['result_sha256'], paper.digest(out))

    def test_reference_files_and_checkpoint_match_paper(self):
        self.assertTrue(all(x['status'] == 'ok' for x in paper.inventory('reference')))
        cfg = json.loads((ROOT / 'configs/paper_assets.json').read_text())
        self.assertEqual(cfg['inputs']['route_b_model']['sha256'],
                         '112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7')
        cfg = json.loads((ROOT / 'configs/paper_tasks.json').read_text())
        self.assertEqual(cfg['method']['members'], ['direct_pc1', 'aligned_c_ensemble', 'b2_verbatim'])
        self.assertEqual(cfg['method']['rrf_k'], 20)
        comparison = json.loads((ROOT / 'configs/v7_b2_joint_formal.json').read_text())
        self.assertNotIn('reference_axis_esm2_mean', comparison['trackA'])

    def test_changed_reference_is_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference = tmp / 'data/paper'; reference.mkdir(parents=True)
            (reference / 'results.csv').write_bytes(b'changed')
            (reference / 'assets.json').write_text(json.dumps({'reference': [{'path': 'data/paper/results.csv',
                'size': 8, 'sha256': hashlib.sha256(b'expected').hexdigest()}]}))
            with patch.object(paper, 'ROOT', tmp), patch.object(paper, 'REFERENCE', reference):
                with self.assertRaisesRegex(ValueError, 'hash-mismatch'):
                    paper.export(tmp / 'output.csv')
            self.assertFalse((tmp / 'output.csv').exists())

    def test_finite_null_probabilities_include_ties(self):
        import numpy as np
        p = paper.screenwide_p(np.array([0., 2., 4.]), np.zeros(3), np.ones(3), np.array([1., 2., 2., 3.]))
        np.testing.assert_array_equal(p, [1., .8, .2])
        # Opposite directions remain separate; passing one never promotes the other.
        other = paper.screenwide_p(np.array([-4.]), np.zeros(1), np.ones(1), np.array([1., 2., 2., 3.]))
        self.assertEqual(other[0], 1.)

    def test_invalid_nulls_and_calibration_fail(self):
        import numpy as np
        for null in [[], [float('nan')]]:
            with self.assertRaises(ValueError): paper.screenwide_p([1.], [0.], [1.], np.array(null))
        for scale in [0., -1., float('nan')]:
            with self.assertRaises(ValueError): paper.screenwide_p([1.], [0.], [scale], [1., 2.])

    def test_output_reuse_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'sentinel.txt').write_text('preserve')
            with self.assertRaisesRegex(ValueError, 'new or empty'): paper.new_output(tmp)
            self.assertEqual((Path(tmp) / 'sentinel.txt').read_text(), 'preserve')

    def test_installer_validates_every_member_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); archive = tmp / 'input.tar.gz'; root = tmp / 'checkout'; root.mkdir()
            expected = b'valid'
            entries = [{'path': name, 'size': len(expected), 'sha256': hashlib.sha256(expected).hexdigest()}
                       for name in ['raw/paper/one', 'raw/paper/two']]
            with tarfile.open(archive, 'w:gz') as t:
                for name, payload in [('raw/paper/one', expected), ('raw/paper/two', b'WRONG')]:
                    info = tarfile.TarInfo(name); info.size = len(payload); t.addfile(info, io.BytesIO(payload))
            with patch.object(paper, 'ROOT', root), patch.object(paper, 'manifest', return_value={'inference': entries}):
                with self.assertRaisesRegex(ValueError, 'checksum mismatch'): paper.install(archive)
            self.assertFalse((root / 'raw/paper/one').exists())

    def test_installer_rejects_paths_outside_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / 'input.tar.gz'
            with tarfile.open(archive, 'w:gz') as t:
                info = tarfile.TarInfo('../escape'); info.size = 1; t.addfile(info, io.BytesIO(b'x'))
            with patch.object(paper, 'manifest', return_value={'inference': []}):
                with self.assertRaisesRegex(ValueError, 'unregistered'): paper.install(archive)
            self.assertFalse((Path(tmp) / 'escape').exists())

    def test_installer_preserves_different_local_asset(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); root = tmp / 'checkout'; root.mkdir(); (root / 'asset').write_bytes(b'local')
            archive = tmp / 'input.tar.gz'; value = b'expected'
            with tarfile.open(archive, 'w:gz') as t:
                info = tarfile.TarInfo('asset'); info.size = len(value); t.addfile(info, io.BytesIO(value))
            entry = {'path': 'asset', 'size': len(value), 'sha256': hashlib.sha256(value).hexdigest()}
            with patch.object(paper, 'ROOT', root), patch.object(paper, 'manifest', return_value={'inference': [entry]}):
                with self.assertRaisesRegex(ValueError, 'different local asset'): paper.install(archive)
            self.assertEqual((root / 'asset').read_bytes(), b'local')


if __name__ == '__main__':
    unittest.main()
