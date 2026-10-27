#!/usr/bin/env python3
"""Compare a completed recipe replay with the published model tensors and gene table."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import paper_release as paper


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True,help='Completed train.py output directory')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    import numpy as np
    import torch
    paper.require('inference')
    run=json.loads((a.run/'run.json').read_text())
    if run['status']!='complete' or run['mode']!='training':p.error('Run is not completed full training')
    if run['implementation']!='historical-v1':p.error('Paper comparison requires historical-v1')
    out=paper.new_output(a.output);reference=ROOT/'raw/tier1_models/route_b';fresh=a.run/'model'
    # Both files have explicit local provenance; reference integrity is checked above.
    x=torch.load(reference/'final_model.pt',map_location='cpu',weights_only=False)['state_dict']
    y=torch.load(fresh/'final_model.pt',map_location='cpu',weights_only=False)['state_dict']
    missing=sorted(set(x)-set(y));extra=sorted(set(y)-set(x));checks=[]
    for key in sorted(set(x)&set(y)):
        u,v=x[key],y[key]
        row=dict(tensor=key,shape_match=list(u.shape)==list(v.shape),dtype_match=str(u.dtype)==str(v.dtype))
        row['exact']=bool(row['shape_match'] and row['dtype_match'] and torch.equal(u,v))
        if row['shape_match']:
            delta=(u.to(torch.float64)-v.to(torch.float64)).abs()
            row.update(max_absolute_error=float(delta.max()) if delta.numel() else 0.,mean_absolute_error=float(delta.mean()) if delta.numel() else 0.)
        checks.append(row)
    old=np.load(reference/'gene_table_final.npy').astype(np.float64);new=np.load(fresh/'gene_table_final.npy').astype(np.float64)
    if old.shape!=new.shape:p.error('Learned gene-table shape differs')
    cosine=np.sum(old*new,axis=1)/(np.linalg.norm(old,axis=1)*np.linalg.norm(new,axis=1))
    table=dict(shape=list(old.shape),exact=bool(np.array_equal(old,new)),max_absolute_error=float(np.max(np.abs(old-new))),mean_absolute_error=float(np.mean(np.abs(old-new))),mean_row_cosine=float(np.mean(cosine)),min_row_cosine=float(np.min(cosine)),reference_sha256=paper.digest(reference/'gene_table_final.npy'),replay_sha256=paper.digest(fresh/'gene_table_final.npy'))
    exact=not missing and not extra and all(r['exact'] for r in checks) and table['exact']
    report=dict(status='exact' if exact else 'different',full_training_completed=True,epochs=len(json.loads((fresh/'epoch_metrics.json').read_text())),missing_tensors=missing,extra_tensors=extra,n_tensors=len(checks),n_exact_tensors=sum(r['exact'] for r in checks),gene_table=table,tensors=checks,training_record_sha256=paper.digest(a.run/'run.json'),script_sha256=paper.digest(Path(__file__)))
    paper.write_json(out/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='tensors'},indent=2))
    return 0 if exact else 2

if __name__=='__main__':raise SystemExit(main())
