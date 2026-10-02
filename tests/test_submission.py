"""Behavioral checks for candidate provenance, portability, and pipeline integration."""

import hashlib
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SubmissionTests(unittest.TestCase):
    def test_published_export_is_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "results.csv"
            subprocess.run(
                [sys.executable, str(ROOT / "predict.py"), "--output", str(out)],
                check=True,
                cwd=tmp,
            )
            self.assertEqual(out.read_bytes(), (ROOT / "results.csv").read_bytes())
            meta = json.loads(out.with_suffix(".metadata.json").read_text())
            self.assertEqual(meta["rows"], 80)
            self.assertEqual(meta["mode"], "published-export")
            self.assertTrue(meta["sources"])
            with out.open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual({r["剂量单位"] for r in rows}, {"micromolar"})

    def test_demo_runs_without_any_historical_results_or_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "clean"
            root.mkdir()
            for name in ["main.py", "predict.py"]:
                shutil.copyfile(ROOT / name, root / name)
            shutil.copytree(
                ROOT / "scripts",
                root / "scripts",
                ignore=shutil.ignore_patterns("__pycache__"),
            )
            shutil.copytree(ROOT / "data/demo", root / "data/demo")
            out = root / "output"
            subprocess.run(
                [
                    sys.executable,
                    str(root / "main.py"),
                    "--mode",
                    "demo",
                    "--output",
                    str(out),
                ],
                check=True,
                cwd=tmp,
            )
            report = json.loads((out / "run.json").read_text())
            self.assertEqual(report["status"], "complete")
            self.assertFalse((root / "results_a6000").exists())
            self.assertFalse((root / "raw").exists())
            metadata = json.loads((out / "results.metadata.json").read_text())
            self.assertEqual(metadata["mode"], "synthetic-demo")
            text = (out / "results.csv").read_text()
            self.assertIn("DEMO-", text)
            self.assertIn("合成演示，非科研候选", text)
            # A planted synthetic positive control must exercise real residual computation.
            self.assertTrue(list((out / "family/tasks_resid").glob("yeast_task_*.tsv")))
            self.assertTrue(
                any(
                    x["name"].startswith("resid_") and "log_sha256" in x
                    for x in report["steps"]
                )
            )
            # Reusing an output directory must fail before stale results can be consumed.
            result = subprocess.run(
                [
                    sys.executable,
                    str(root / "main.py"),
                    "--mode",
                    "demo",
                    "--output",
                    str(out),
                ],
                capture_output=True,
                text=True,
                cwd=tmp,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("new or empty", result.stderr)

    def test_asset_installer_rejects_wrong_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            wrong = Path(tmp) / "model.pt"
            wrong.write_bytes(b"not the frozen model")
            proc = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts/install_assets.py"),
                    "--weights",
                    str(wrong),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("checksum mismatch", proc.stderr)

    def test_frozen_row_order_bytes(self):
        p = ROOT / "raw/tier1_models/route_b/A2_init.tsv"
        if not p.exists():
            self.skipTest("large assets omitted from review package")
        self.assertEqual(
            hashlib.sha256(p.read_bytes()).hexdigest(),
            "411e09020354463f986a47f8f8f48dbae0ecd547df8719e7c2032ad37a6e5c7c",
        )


if __name__ == "__main__":
    unittest.main()


class PreprocessingTests(unittest.TestCase):
    def test_response_header_variants_and_output_guard(self):
        sys.path.insert(0, str(ROOT))
        from scripts.preprocess_hiphop import build_hiphop_response, HiphopResponseError
        import numpy as np

        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            with (base / "sdrf.tsv").open("w") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(
                    [
                        "Array Data File",
                        "Comment[inchikey]",
                        "Factor Value[dose]",
                        "Unit[concentration unit]",
                    ]
                )
                w.writerow(["01_vehicle", "", "1", "percent"])
                w.writerow(["01_drug1", "TEST1", "2", "micromolar"])
                w.writerow(["01_drug2", "TEST2", "3", "nanomolar"])
            outputs = []
            for explicit_header in [False, True]:
                expr = base / f"expr{explicit_header}.gz"
                with gzip.open(expr, "wt") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow(
                        ([""] if explicit_header else [])
                        + ["01_vehicle", "01_drug1", "01_drug2"]
                    )
                    w.writerows(
                        [["Y1.a", 10, 12, 14], ["Y1.b", 20, 24, 26], ["Y2.a", 8, 5, 10]]
                    )
                name = f"out{explicit_header}.npz"
                record = f"out{explicit_header}.json"
                build_hiphop_response(
                    tmp,
                    expr_path=str(expr),
                    sdrf_path="sdrf.tsv",
                    output_npz=name,
                    output_json=record,
                )
                with np.load(base / name) as z:
                    np.testing.assert_array_equal(z["z_score"], [[0, 3, 5], [0, -3, 2]])
                    self.assertEqual(
                        z["dose_units"].tolist(), ["percent", "micromolar", "nanomolar"]
                    )
                before = (base / name).read_bytes()
                with self.assertRaises(HiphopResponseError):
                    build_hiphop_response(
                        tmp,
                        expr_path=str(expr),
                        sdrf_path="sdrf.tsv",
                        output_npz=name,
                        output_json=record,
                    )
                self.assertEqual(before, (base / name).read_bytes())
