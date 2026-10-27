#!/usr/bin/env python3
"""Run and verify the manuscript's functional, backbone and representation comparisons."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['core','representation-grid','routes'],default='core')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args()
    if a.threads<1:p.error('--threads must be positive')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]=str(a.threads)
    import pandas as pd
    from scripts import paper_release as paper
    paper.require('inference');paper.require('comparisons')
    if a.mode=='routes':
        return subprocess.call([sys.executable,str(ROOT/'scripts/paper_route_comparisons.py'),'--output',str(a.output),'--threads',str(a.threads)])
    output=paper.new_output(a.output);logs=paper.new_output(output/'logs')
    commands=[]
    if a.mode=='core':
        for name,script,config in [('functional_v9','v9_esm2_basic_run.py','v9_esm2_basic.json'),('b2_comparison','v7_b2_formal_run.py','v7_b2_joint_formal.json'),('go_comparison','v5_joint_formal_run.py','v5_1_joint_formal.json')]:
            commands.append((name,[sys.executable,str(ROOT/'scripts'/script),'--config',str(ROOT/'configs'/config),'--output',str(output/name)],{}))
        commands.extend([
            ('functional_ev1a',[sys.executable,str(ROOT/'scripts/ev1a_run.py')],{'YB_EV1A_OUT':str(output/'functional_ev1a')}),
            ('backbone',[sys.executable,str(ROOT/'scripts/model_selection/step1_scF_selection/norman_retrieval.py')],{'YB_MODEL_SELECTION_OUTPUT':str(output/'backbone')})])
    else:
        b=ROOT/'raw/model_selection/representation_grid'
        commands=[('representation_grid',[sys.executable,str(ROOT/'scripts/v5_joint_formal_run.py'),'--config',str(ROOT/'configs/v5_1_joint_formal.json'),'--repr-grid',str(output/'representation_grid'),'--scf-npy',str(b/'scfoundation.npy'),'--scyeast-prior-npy',str(b/'scyeast_prior.npy'),'--raw-h5ad',str(b/'observation_index.tsv')],{})]
    for name,command,extra in commands:
        print('Running',name,flush=True)
        with (logs/(name+'.log')).open('w') as log:
            subprocess.run(command,cwd=ROOT,env=dict(os.environ,**extra),stdout=log,stderr=subprocess.STDOUT,check=True)
    reference=ROOT/'raw/model_selection/reference';checks=[]
    if a.mode=='core':
        pairs=[('functional_v9/per_query_rr_basic_structure.tsv','v9_rr.tsv'),('functional_ev1a/per_query_rr_ev1a.tsv','ev1a_rr.tsv'),('b2_comparison/per_query_ranks_trackA_b2_primary.tsv','b2_ranks.tsv'),('go_comparison/per_query_ranks_trackA.tsv','go_ranks.tsv')]
        for name,ref in pairs:
            x=pd.read_csv(output/name,sep='\t');y=pd.read_csv(reference/'core'/ref,sep='\t')
            pd.testing.assert_frame_equal(x,y,check_exact=True)
            checks.append(dict(output=name,rows=len(x),columns=len(x.columns),values_exact=True,reference_sha256=paper.digest(reference/'core'/ref)))
        x=json.loads((output/'backbone/retrieval_result.json').read_text());y=json.loads((reference/'step1_scF_selection/retrieval_result.json').read_text())
        if x!=y:raise ValueError('Backbone comparison differs')
        pooled=pd.concat([pd.read_csv(output/name,sep='\t') for name,_ in pairs[:2]])
        mrr={c:float(pooled[c].mean()) for c in ['rr_human_only','rr_yeast_only','rr_comparator_best_single','rr_joint_nested']}
        report=dict(status='passed',checks=checks,backbone_all_fields_equal=True,pooled_n=len(pooled),pooled_mrr=mrr)
    else:
        x=json.loads((output/'representation_grid/repr_grid_result.json').read_text());y=json.loads((reference/'step2_scyeast_exclusion/repr_grid_result.json').read_text())
        if x['cells']!=y['cells']:raise ValueError('Representation-grid comparison differs')
        report=dict(status='passed',all_grid_cells_equal=True,cells=x['cells'])
    report['script_sha256']=paper.digest(Path(__file__))
    paper.write_json(output/'verification.json',report);print(json.dumps(report,indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
