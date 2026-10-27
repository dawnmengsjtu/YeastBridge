#!/usr/bin/env python3
"""Rerun the paper's frozen A2/B2/C2/scYeast T2–T5 comparisons."""
import argparse
import json
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args()
    if a.threads < 1:p.error('--threads must be positive')
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]=str(a.threads)
    import numpy as np
    import pandas as pd
    from scripts.paper_release import new_output,write_json,digest,require
    require('comparisons')
    from scripts.model_selection.step4_route_confirmation.eval import tasks
    output=new_output(a.output)
    reference=ROOT/'raw/model_selection/reference/step4_route_confirmation'
    ref=json.loads((reference/'five_route_results.json').read_text())
    features=['A2_scf_ortholog','B2_scf_esm2inject','C2_scf_scratch','routed_scyeast']
    routes=ROOT/'raw/model_selection/routes/transfer_routes/scf_routes'
    for route,feature in zip(['A2','B2','C2'],features):
        if digest(routes/route/'gene_table_final.npy') != ref['asset_hashes'][feature]:
            raise ValueError('Route table differs from the registered comparison: '+feature)
    results={};differences=[]
    for feature in features:
        results[feature]={}
        for task,fn in [('T2',tasks.run_t2),('T3',tasks.run_t3),('T4',tasks.run_t4),('T5',tasks.run_t5)]:
            metrics,detail=fn(feature=feature,seed=42)
            results[feature][task]=metrics
            detail.to_csv(output/f'{task.lower()}_{feature}.csv',index=False)
            for k,v in ref['tasks'][feature][task].items():
                if isinstance(v,(int,float)) and abs(metrics[k]-v)>1e-6:
                    differences.append(dict(feature=feature,task=task,metric=k,expected=v,actual=metrics[k]))
            print(feature,task,json.dumps(metrics),flush=True)
    detail_checks={}
    for task in ['t3','t5']:
        fresh=pd.read_csv(output/f'{task}_B2_scf_esm2inject.csv')
        frozen=pd.read_csv(reference/f'{task}_B2_scf_esm2inject_ridge_s42_detail.csv')
        pd.testing.assert_frame_equal(fresh,frozen,check_exact=False,atol=1e-6,rtol=0)
        detail_checks[task]=dict(rows=len(fresh),all_fields_equal_within=1e-6,reference_sha256=digest(reference/f'{task}_B2_scf_esm2inject_ridge_s42_detail.csv'))
    write_json(output/'result.json',results)
    write_json(output/'verification.json',dict(status='passed' if not differences else 'different',
               differences=differences,B2_detail_checks=detail_checks,scope='Frozen representations; T2–T5 evaluation retrained from fixed data splits',
               excluded_independent_direct_ESM2=True,script_sha256=digest(Path(__file__))))
    return bool(differences)

if __name__=='__main__':raise SystemExit(main())
