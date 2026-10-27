#!/usr/bin/env python3
"""Fresh candidate/native docking from the manuscript's prepared coordinates.

Source coordinates are read-only. Search outputs and caches are never copied.
The complete control-panel and contact interpretation are separate analyses.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paper_release import digest, new_output, write_json
GNINA_SHA256 = '5d33538324b40050a03aa262d51832837e0ea6cc100945abbd2d7b732589690e'


def dock(job):
    case, original = job
    os.environ['YB_DOCK_CASE_DIR'] = case
    path = ROOT / 'scripts/structural/dock_run_eval.py'
    spec = importlib.util.spec_from_file_location('docking_case', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return dict(case=case, **module.ccr4_dock_one(tuple(original)))


def compare_searches(rows, source):
    """Report every search, including missing references and nonfinite scores."""
    import pandas as pd
    comparisons, missing, tables = [], [], {}
    for row in rows:
        relative = row['case'].split('/cases/')[-1]
        identity = {key: row[key] for key in ['structure', 'ligand', 'seed', 'native', 'status']}
        identity['case'] = relative
        unit=row.get('score_unit')
        identity['score_unit']=str(unit) if pd.notna(unit) else 'kcal/mol'
        ref = source / relative / ('native_docking.tsv' if row['native'] else 'panel_docking.tsv')
        if not ref.exists():
            missing.append(dict(**identity, reason='Reference table absent'))
            continue
        if ref not in tables:
            tables[ref] = pd.read_csv(ref, sep='\t')
        table = tables[ref]
        match = table[table.structure.eq(row['structure']) & table.ligand.eq(row['ligand']) & table.seed.eq(row['seed'])]
        if len(match) != 1:
            missing.append(dict(**identity, reason='Reference row missing or ambiguous', matches=len(match)))
            continue
        original = match.iloc[0]
        errors, nonfinite = {}, []
        for key in ['affinity', 'cnn_affinity', 'cnn_score', 'top_pose_rmsd']:
            if key not in row or key not in original:
                continue
            x, y = float(row[key]), float(original[key])
            if math.isfinite(x) and math.isfinite(y):
                errors[key] = abs(x - y)
            elif not ((math.isnan(x) and math.isnan(y)) or x == y):
                nonfinite.append(key)
        comparisons.append(dict(**identity, absolute_errors=errors, nonfinite_mismatches=nonfinite))
    columns = sorted({(key,row['score_unit'] if key=='affinity' else None)
                      for row in comparisons for key in row['absolute_errors']})
    summary = {}
    for key,unit in columns:
        values = [row['absolute_errors'][key] for row in comparisons if key in row['absolute_errors']
                  and (unit is None or row['score_unit']==unit)]
        label=key+' ['+unit+']' if unit else key
        summary[label] = dict(compared=len(values), exact=sum(v == 0 for v in values),
                            max_absolute_error=max(values), mean_absolute_error=sum(values)/len(values))
    qualification = []
    native_rows = pd.DataFrame([row for row in rows if row['native']])
    if len(native_rows):
        for (case, sid), group in native_rows.groupby(['case', 'structure']):
            relative = case.split('/cases/')[-1]
            ref = source / relative / 'redocking_qualification.tsv'
            if not ref.exists():
                continue
            table = pd.read_csv(ref, sep='\t')
            original = table[table.structure.eq(sid)]
            if len(original) != 1:
                continue
            original = original.iloc[0]
            values = group.loc[group.status.eq('ok'), 'top_pose_rmsd']
            median = float(values.median())
            plan = json.loads((source / relative / 'analysis_plan.json').read_text())
            qualified = len(values) == len(plan['seeds']) and median <= 2.0
            qualification.append(dict(case=relative, structure=sid, n_success=len(values),
                median_rmsd=median, reference_median_rmsd=float(original.median_top_pose_rmsd),
                qualified=qualified, reference_qualified=bool(original.qualified),
                qualification_match=qualified == bool(original.qualified)))
    return dict(comparisons=comparisons, reference_rows_compared=len(comparisons),
                missing_references=missing, numerical_summary=summary,
                native_qualification_comparisons=qualification,
                nonfinite_mismatch_jobs=sum(bool(row['nonfinite_mismatches']) for row in comparisons))


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='structural_current from Supplementary Data 1')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--gnina', type=Path)
    parser.add_argument('--gnina-launcher-json', default='[]')
    parser.add_argument('--gnina-devices', default='0')
    parser.add_argument('--plan-only', action='store_true', help='Prepare and verify inputs without running searches')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if not args.gnina or not args.gnina.is_file():
        parser.error('--gnina must name the GNINA 1.3.2 executable')
    if digest(args.gnina) != GNINA_SHA256:
        parser.error('GNINA checksum differs; use the official v1.3.2 asset gnina.1.3.2')
    launcher = json.loads(args.gnina_launcher_json)
    if not isinstance(launcher,list) or not all(isinstance(x,str) for x in launcher):
        parser.error('launcher must be a JSON array of command arguments')
    source = args.source.resolve()
    output = new_output(args.output)
    os.environ['YB_GNINA_DEVICES'] = args.gnina_devices
    verdict = pd.read_csv(source/'pair_verdicts.tsv',sep='\t',keep_default_na=False)
    cases = sorted({r.evaluation_case.split('/structural_validation/')[-1] for r in verdict.itertuples() if r.docking_tested})
    cases.append('CLIC1/gnina_refinement/search_depth')
    jobs, manifests, skipped = [], [], []
    for relative in cases:
        original = source/relative
        case = output/'cases'/relative
        case.mkdir(parents=True)
        plan = json.loads((original/'analysis_plan.json').read_text())
        plan['gnina_binary'] = str(args.gnina.resolve())
        plan['gnina_launcher'] = launcher
        write_json(case/'analysis_plan.json',plan)
        for name in ['structures.json','ligand_preparation.tsv','stereochemistry_feasibility.tsv']:
            if (original/name).is_file():
                shutil.copyfile(original/name,case/name)
                manifests.append(dict(path=str((original/name).relative_to(source)),sha256=digest(original/name)))
        for folder in ['receptors','structures','ligands']:
            (case/folder).mkdir()
            for p in (original/folder).glob('*'):
                if p.is_file() and p.suffix in ['.pdb','.pdbqt','.sdf']:
                    shutil.copyfile(p,case/folder/p.name)
                    manifests.append(dict(path=str(p.relative_to(source)),sha256=digest(p)))
        (case/'poses').mkdir()
        ligands = pd.read_csv(case/'ligand_preparation.tsv',sep='\t')
        roles = ['candidate','experimental_control'] if relative.endswith('search_depth') else ['candidate']
        ligands = ligands[ligands.role.isin(roles)]
        if (case/'stereochemistry_feasibility.tsv').exists():
            feasible = pd.read_csv(case/'stereochemistry_feasibility.tsv',sep='\t').set_index('ligand').geometry_feasible
            ligands = ligands[ligands.ligand.map(feasible).fillna(False)]
        complete = ligands.groupby('inchikey').status.apply(lambda x:x.eq('ok').all())
        ligands = ligands[ligands.inchikey.isin(complete[complete].index)]
        boxes = json.loads((case/'structures.json').read_text())
        for sid in plan['structures']:
            box = next(b for b in boxes if b['structure']==sid)
            if box.get('native_available',True) and (case/'ligands'/f'{sid}_native.pdbqt').exists():
                jobs.extend((str(case),(sid,sid+'_native',seed,True)) for seed in plan['seeds'])
            else:
                skipped.append(dict(case=relative,structure=sid,operation='native',reason='No prepared native ligand'))
            jobs.extend((str(case),(sid,r.ligand,seed,False)) for r in ligands.itertuples() for seed in plan['seeds'])
    plan = dict(cases=cases,jobs=len(jobs),source_manifest=manifests,skipped=skipped,
                source_poses_reused=False,source_docking_cache_reused=False,
                preparation='Frozen receptor and ligand coordinates',
                gnina_binary_sha256=GNINA_SHA256,
                script_sha256=digest(Path(__file__)),engine_script_sha256=digest(ROOT/'scripts/structural/dock_run_eval.py'))
    plan['jobs_list']=[dict(case=case.split('/cases/')[-1],structure=job[0],ligand=job[1],seed=job[2],native=job[3]) for case,job in jobs]
    for case,job in jobs:
        sid,ligand,seed,native=job
        required=[Path(case)/'receptors'/(sid+'.pdbqt'),Path(case)/'ligands'/(ligand+'.pdbqt'),Path(case)/'ligands'/(ligand+'.sdf')]
        if native:required.append(Path(case)/'ligands'/(ligand+'_reference.sdf'))
        for path in required:
            if not path.is_file():raise ValueError('Missing prepared input: '+str(path))
    write_json(output/'execution_plan.json',plan)
    if args.plan_only:
        print(json.dumps(dict(status='prepared',jobs=len(jobs),searches_executed=0)))
        return 0
    started=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(dock,job) for job in jobs]):
            row=future.result();rows.append(row)
            pd.DataFrame(rows).to_csv(output/'fresh_docking.tsv',sep='\t',index=False)
            print(json.dumps(dict(completed=len(rows),total=len(jobs),case=row['case'].split('/cases/')[-1],status=row['status'])),flush=True)
    comparison = compare_searches(rows, source)
    failed=[r for r in rows if r['status']!='ok']
    write_json(output/'verification.json',dict(status='complete' if not failed else 'failed',fresh_searches=len(rows),
        successful=len(rows)-len(failed),failed=failed,**comparison,seconds=time.monotonic()-started,
        scope='Prepared-coordinate candidate/native redocking; complete background/control-panel reruns are separate'))
    return bool(failed)

if __name__=='__main__':
    raise SystemExit(main())
