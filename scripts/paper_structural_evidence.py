#!/usr/bin/env python3
"""Recompute contact evidence and manuscript prioritization from archived poses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.paper_release import digest,new_output,write_json
ARCHIVE_SHA='486c3c896ed7709c9eda3b075c3d7d039290fa5a5a38550a866631d33755da1a'
ORIGINAL='/public/home/mengxl/dzy/YeastBridge_0912/results_a6000/ion_channel_scope_20261012/full_measured_assay/structural_validation'


def rebase(value,root):
    if isinstance(value,str):return value.replace(ORIGINAL,str(root))
    if isinstance(value,list):return [rebase(v,root) for v in value]
    if isinstance(value,dict):return {k:rebase(v,root) for k,v in value.items()}
    return value


def evidence(root):
    import pandas as pd
    control=root/'control_calibration'
    comparison=json.loads((control/'reference_pocket_comparison.json').read_text())
    clic=json.loads((root/'CLIC1/gnina_refinement/search_depth/literature_contact_evidence.json').read_text())
    models=json.loads((control/'interaction_calibration.json').read_text())
    positive={}
    for row in comparison['records']:
        if row['candidate_id'] and row['target']=='RYR1':
            positive[row['candidate_id']]=all(row[k].get('TRP996.A',0)==row['n_seeds'] for k in ['hotspot_contact_counts','cnn_hotspot_contact_counts'])
    comparator=next(r for r in clic['records'] if r['role']!='candidate' and r['inchikey']!=clic['native_inchikey'])
    for row in clic['records']:
        if row['role']=='candidate':
            positive[row['name']]=all(all(s['hotspot_contacts'].get(k) for k in ['TYR214.A','PHE83.B']) for s in row['seeds']) and row['median_cnn_affinity']>comparator['median_cnn_affinity']
    primary=next(m for m in models['models'] if m['method']=='contact_only')
    for row in primary['candidates']:
        positive[row['candidate_id']]=min([row['score']]+[s['score'] for s in row['seed_scores']]+[s['score'] for s in row['study_omission_scores']])>=.5
    rows=[]
    for row in pd.read_csv(root/'pair_verdicts.tsv',sep='\t').itertuples():
        good=bool(row.protocol_qualified and row.n_geometry_pass>=2 and positive.get(row.candidate_id,False))
        rows.append(dict(candidate_id=row.candidate_id,target_id=row.target_id,
                         docking_status='Positive' if good else 'Not positive' if row.docking_tested else 'Not tested'))
    return pd.DataFrame(rows)


def scientific_payload(value):
    # Runtime paths, input file hashes and library-version records differ after rebasing.
    if isinstance(value,dict):
        value=dict(value)
        # Older Vina contact records left their common kcal/mol unit implicit.
        if 'seed' in value and 'size_corrected_affinity' in value and 'score_unit' not in value:
            value['score_unit']='kcal/mol'
        return {k:scientific_payload(v) for k,v in value.items() if k not in ['sources','source_sha256','case','path','specification','sklearn_version','prolif_version','software']}
    if isinstance(value,list):return [scientific_payload(v) for v in value]
    return value


def install_fresh_searches(root, replay):
    """Overlay verified fresh candidate/native searches; preserve background controls."""
    import pandas as pd
    verification=json.loads((replay/'verification.json').read_text())
    plan=json.loads((replay/'execution_plan.json').read_text())
    if verification['status']!='complete' or verification['failed']:
        raise ValueError('Fresh docking has not completed successfully')
    rows=pd.read_csv(replay/'fresh_docking.tsv',sep='\t')
    if len(rows)!=plan['jobs'] or not rows.status.eq('ok').all():
        raise ValueError('Fresh docking row coverage differs')
    for original,name in [('verification.json','redocking_verification_input.json'),
                          ('execution_plan.json','redocking_execution_plan_input.json'),
                          ('fresh_docking.tsv','redocking_table_input.tsv')]:
        shutil.copyfile(replay/original,root.parent/name)
    required=set()
    for row in rows.itertuples():
        relative=row.case.split('/cases/')[-1]
        required.update(relative+'/'+name for name in ['structures.json','ligand_preparation.tsv',
            'receptors/'+row.structure+'.pdbqt','ligands/'+row.ligand+'.sdf','ligands/'+row.ligand+'.pdbqt'])
        if row.native:required.add(relative+'/ligands/'+row.ligand+'_reference.sdf')
    verified=set();unused_absent=0
    for row in plan['source_manifest']:
        path=root/row['path']
        if not path.exists() and row['path'] not in required:
            unused_absent+=1
            continue
        if not path.exists() or digest(path)!=row['sha256']:
            raise ValueError('Prepared docking input differs: '+row['path'])
        verified.add(row['path'])
    if not required.issubset(verified):raise ValueError('A required docking input lacks verified provenance')
    for case,group in rows.groupby('case'):
        relative=case.split('/cases/')[-1]
        if relative not in plan['cases']:raise ValueError('Unexpected docking case')
        target=root/relative;fresh=replay/'cases'/relative
        original_plan=json.loads((target/'analysis_plan.json').read_text())
        executed_plan=json.loads((fresh/'analysis_plan.json').read_text())
        runtime_keys={'gnina_binary','gnina_launcher'}
        if {k:v for k,v in original_plan.items() if k not in runtime_keys} != {k:v for k,v in executed_plan.items() if k not in runtime_keys}:
            raise ValueError('Docking protocol changed: '+relative)
        for file in (fresh/'poses').iterdir():
            if file.is_file():shutil.copyfile(file,target/'poses'/file.name)
        for native,name in [(False,'panel_docking.tsv'),(True,'native_docking.tsv')]:
            replacement=group.loc[group.native.eq(native)].drop(columns='case').dropna(axis=1,how='all')
            if replacement.empty:continue
            path=target/name
            if path.exists():
                old=pd.read_csv(path,sep='\t')
                keys=['structure','ligand','seed']
                selected=pd.MultiIndex.from_frame(replacement[keys])
                old=old.loc[~pd.MultiIndex.from_frame(old[keys]).isin(selected)]
                replacement=pd.concat([old,replacement],ignore_index=True)
            replacement.to_csv(path,sep='\t',index=False)
        native_rows=group[group.native]
        if len(native_rows):
            case_plan=json.loads((target/'analysis_plan.json').read_text());checks=[]
            for sid,native in native_rows.groupby('structure'):
                median=float(native.top_pose_rmsd.median())
                checks.append(dict(structure=sid,n_success=len(native),median_top_pose_rmsd=median,
                    reference_type=case_plan.get('reference_type','experimental_native_pose'),
                    qualified=len(native)==len(case_plan['seeds']) and median<=2.0))
            prior=pd.read_csv(target/'redocking_qualification.tsv',sep='\t').set_index('structure')
            if any(row['structure'] not in prior.index or
                   row['qualified'] != bool(prior.loc[row['structure'],'qualified']) for row in checks):
                raise ValueError('Native qualification changed; the full pair-verdict analysis must be rerun before prioritization')
            pd.DataFrame(checks).to_csv(target/'redocking_qualification.tsv',sep='\t',index=False)
    return dict(fresh_candidate_native_searches=len(rows),background_control_searches='frozen',
                required_prepared_inputs_verified=len(required),unused_background_input_files_absent=unused_absent,
                replay_verification_sha256=digest(root.parent/'redocking_verification_input.json'),
                replay_table_sha256=digest(root.parent/'redocking_table_input.tsv'))


def main():
    import pandas as pd
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--redocking',type=Path,help='Completed paper_redocking.py output to audit new candidate/native poses')
    a=p.parse_args()
    if digest(a.archive)!=ARCHIVE_SHA:p.error('Source archive is not the frozen manuscript release')
    output=new_output(a.output);root=output/'structural_current';references={}
    stages=[('control_calibration','reference_pocket_contacts','reference_pocket_contacts.json'),
            ('control_calibration','reference_pocket_cnn_contacts','reference_pocket_cnn_contacts.json'),
            ('control_calibration','reference_pocket_summary','reference_pocket_comparison.json'),
            ('CLIC1/gnina_refinement/search_depth','literature_contacts','literature_contact_evidence.json'),
            ('control_calibration','interaction_fingerprints','interaction_fingerprints.json'),
            ('control_calibration','interaction_calibration','interaction_calibration.json')]
    with zipfile.ZipFile(a.archive) as archive:
        for item in archive.infolist():
            if not item.filename.startswith('structural_current/') or item.is_dir():continue
            path=output/item.filename
            if not path.resolve().is_relative_to(root.resolve()):raise ValueError('Unsafe archive path')
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(archive.read(item))
        expected=pd.read_csv(archive.open('manuscript_data/docking_current_evidence.csv'))
    for case,stage,name in stages:
        references[f'{case}/{name}']=json.loads((root/case/name).read_text())
    replay=install_fresh_searches(root,a.redocking.resolve()) if a.redocking else None
    # Only this private, new work directory is rebased or updated.
    for path in root.rglob('*.json'):
        value=json.loads(path.read_text());write_json(path,rebase(value,root))
    logs=new_output(output/'logs');checks={}
    for case,stage,name in stages:
        print('Recomputing',stage,flush=True)
        with (logs/(stage+'.log')).open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/structural/dock_run_eval.py'),'--mode','ccr4'],
                env=dict(os.environ,YB_DOCK_CASE_DIR=str(root/case),YB_DOCK_STAGE=stage,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4'),
                stdout=log,stderr=subprocess.STDOUT,check=True)
        result=json.loads((root/case/name).read_text())
        key=f'{case}/{name}'
        checks[key]=dict(scientific_payload_equal=scientific_payload(result)==scientific_payload(references[key]),
                        result_sha256=digest(root/case/name))
        write_json(output/'stage_checks.json',checks)
    if not a.redocking and not all(c['scientific_payload_equal'] for c in checks.values()):
        raise ValueError('Recomputed structural evidence differs; inspect stage_checks.json')
    table=evidence(root)
    expected=expected[table.columns]
    statuses_match=table.equals(expected)
    if not a.redocking:pd.testing.assert_frame_equal(table,expected,check_exact=True)
    table.to_csv(output/'docking_current_evidence.csv',index=False)
    selected=table.loc[table.docking_status.eq('Positive'),'candidate_id'].tolist()
    write_json(output/'verification.json',dict(status='passed' if statuses_match else 'different',candidate_statuses_match=statuses_match,
        prioritized=selected,n_prioritized=len(selected),stages=checks,
        fresh_contact_calculation=True,fresh_contact_models=True,fresh_docking=bool(a.redocking),
        docking_replay=replay,
        interpretation='Descriptive pocket-contact/reference prioritization, separate from the stricter binding-control classifier',
        source_archive_sha256=ARCHIVE_SHA,script_sha256=digest(Path(__file__))))
    print(json.dumps(dict(prioritized=selected,checks=checks),indent=2))
    return 0 if statuses_match else 2

if __name__=='__main__':raise SystemExit(main())
