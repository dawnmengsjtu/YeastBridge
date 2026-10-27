"""Scientific invariants for the fresh-null and family-projection entry points."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import paper_release as paper
from scripts import project_io
from scripts.paper_permutations import permute_within_assays
from scripts.paper_calibration import residual_profiles
from scripts.paper_structural_evidence import scientific_payload
from scripts.install_training_assets import install_verified
from scripts.paper_redocking import compare_searches


class FreshReproductionTests(unittest.TestCase):
    def test_source_snapshot_revision_uses_embedded_metadata(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(project_io, 'ROOT', Path(tmp)), patch.object(project_io.subprocess, 'check_output', side_effect=AssertionError('Must not inspect a parent repository')):
            root=Path(tmp)
            self.assertEqual(project_io.revision(), 'source-archive (no embedded code revision)')
            (root/'BUNDLE.json').write_text(json.dumps({'code_revision':'legacy'}))
            self.assertEqual(project_io.revision(), 'legacy')
            (root/'CODE_VERSION.json').write_text(json.dumps({'commit':'a'*40}))
            self.assertEqual(project_io.revision(), 'a'*40)

    def test_comparison_configs_use_declared_release_inputs(self):
        assets={row['path']:row for group in paper.manifest().values() for row in group}
        fields=[('path','sha256'),('path_abs','sha256'),('index_path','index_sha256'),
                ('index_path_abs','index_sha256'),('gene_order_path','gene_order_sha256')]
        for name in ['v5_1_joint_formal.json','v7_b2_joint_formal.json','v9_esm2_basic.json','ev_v7_mirror.json']:
            config=json.loads((ROOT/'configs'/name).read_text())
            for key,entry in config['inputs'].items():
                if not isinstance(entry,dict):continue
                for path_key,sha_key in fields:
                    if path_key not in entry:continue
                    path=entry[path_key]
                    with self.subTest(config=name,input=key,path=path):
                        self.assertIn(path,assets)
                        if entry.get(sha_key):self.assertEqual(entry[sha_key],assets[path]['sha256'])

    def test_redocking_comparison_exposes_missing_and_nonfinite_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp);case=source/'target';case.mkdir()
            pd.DataFrame([dict(structure='PDB',ligand='candidate',seed=42,affinity=-7.0,cnn_score=0.7)]).to_csv(case/'panel_docking.tsv',sep='\t',index=False)
            row=dict(case='/run/cases/target',structure='PDB',ligand='candidate',seed=42,native=False,status='ok',affinity=-6.8,cnn_score=float('nan'))
            report=compare_searches([row,dict(row,ligand='missing')],source)
            self.assertEqual(report['reference_rows_compared'],1)
            self.assertEqual(len(report['missing_references']),1)
            self.assertEqual(report['nonfinite_mismatch_jobs'],1)
            self.assertAlmostEqual(report['numerical_summary']['affinity [kcal/mol]']['max_absolute_error'],0.2)
            self.assertEqual(report['numerical_summary']['affinity [kcal/mol]']['exact'],0)

    def test_structural_comparison_retains_scores_contacts_and_units(self):
        original={'seed':42,'size_corrected_affinity':-7.2,'contacts':['A1'],
                  'source_sha256':'old','path':'old'}
        replay=dict(original,source_sha256='new',path='new',score_unit='kcal/mol')
        self.assertEqual(scientific_payload(original),scientific_payload(replay))
        for key,value in [('score_unit','kJ/mol'),('size_corrected_affinity',-7.1),('contacts',[])]:
            changed=dict(replay,**{key:value})
            self.assertNotEqual(scientific_payload(original),scientific_payload(changed))

    def test_training_install_preserves_different_existing_weights(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';target=root/'target'
            source.write_bytes(b'published model');expected=paper.digest(source)
            target.write_bytes(b'local trained model')
            with self.assertRaisesRegex(ValueError,'Refusing to replace'):
                install_verified(source,target,expected)
            self.assertEqual(target.read_bytes(),b'local trained model')
            fresh=root/'new/model'
            install_verified(source,fresh,expected)
            install_verified(source,fresh,expected)
            self.assertEqual(paper.digest(fresh),expected)

    def test_stratified_permutations_preserve_each_assay(self):
        labels=np.array(['HOP','HIP','HIP','HOP','HOP','HIP','HOP'])
        rng=np.random.default_rng(20261024)
        first=[]
        for _ in range(12):
            p=permute_within_assays(rng,labels);first.append(p)
            np.testing.assert_array_equal(np.sort(p),np.arange(7))
            np.testing.assert_array_equal(labels[p],labels)
        replay=np.random.default_rng(20261024)
        for expected in first:np.testing.assert_array_equal(expected,permute_within_assays(replay,labels))
        self.assertGreater(len({tuple(x) for x in first}),1)

    def family_fixture(self,root,rows,non_estimable):
        profiles=root/'profiles';profiles.mkdir()
        reference=root/'reference';reference.mkdir()
        pd.DataFrame({'target_id':list(rows),'profile_family':['FAM1','FAM1','FAM2']}).to_csv(reference/'family_membership.tsv',sep='\t',index=False)
        (reference/'residual_export_record.json').write_text(json.dumps({'non_estimable_targets':non_estimable}))
        (profiles/'export_record.json').write_text(json.dumps({'tail_strains':[]}))
        for target,values in rows.items():
            pd.DataFrame({'rank':[1,2,3],'yeast_gene':['g1','g2','g3'],'ensemble_score':values}).to_csv(profiles/f'yeast_task_{target}.tsv',sep='\t',index=False)
        return profiles,reference

    def test_residual_projection_removes_family_component(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows={'A':[1,2,3],'B':[-2,1,2],'C':[0,1,2]}
            profiles,reference=self.family_fixture(root,rows,[])
            with patch.object(paper,'REFERENCE',reference):residual_profiles(profiles,root/'out')
            mean=np.mean([rows['A'],rows['B']],axis=0)
            for target in ['A','B']:
                result=pd.read_csv(root/'out'/f'yeast_task_{target}.tsv',sep='\t').set_index('yeast_gene').ensemble_score.reindex(['g1','g2','g3']).to_numpy()
                self.assertLess(abs(result@mean),1e-8)
            self.assertFalse((root/'out/yeast_task_C.tsv').exists())

    def test_identical_family_profiles_remain_non_estimable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);profiles,reference=self.family_fixture(root,{'A':[1,2,3],'B':[1,2,3],'C':[0,1,2]},['A','B'])
            with patch.object(paper,'REFERENCE',reference):self.assertEqual(residual_profiles(profiles,root/'out'),['A','B'])
            for target in ['A','B']:
                values=pd.read_csv(root/'out'/f'yeast_task_{target}.tsv',sep='\t').ensemble_score
                np.testing.assert_array_equal(values,np.zeros(3))

    def test_extended_groups_cannot_replace_core_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'assets.json').write_text(json.dumps({'reference':[]}))
            (root/'extended_assets.json').write_text(json.dumps({'reference':[]}))
            with patch.object(paper,'REFERENCE',root):
                with self.assertRaisesRegex(ValueError,'shadow'):paper.manifest()

if __name__=='__main__':unittest.main()
