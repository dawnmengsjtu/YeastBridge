#!/usr/bin/env python3
"""Rebuild profiles, calibration, maximum nulls and the final manuscript screen."""
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
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--response',type=Path,default=ROOT/'raw/paper/strain_response.npz')
    p.add_argument('--device',default='cuda:0')
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args()
    if a.threads<1:p.error('--threads must be positive')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]=str(a.threads)
    import torch
    from scripts import paper_release as paper
    paper.require('all')
    if a.device.startswith('cuda') and not torch.cuda.is_available():p.error('A CUDA GPU is required for the requested device')
    if paper.digest(a.response)!=paper.digest(ROOT/'raw/paper/strain_response.npz'):p.error('Response checksum differs')
    output=paper.new_output(a.output);logs=paper.new_output(output/'logs')
    stages=[('profiles',[sys.executable,str(ROOT/'main.py'),'--mode','profiles','--threads',str(a.threads),'--output',str(output/'profiles')]),
            ('calibration',[sys.executable,str(ROOT/'scripts/paper_calibration.py'),'--profiles',str(output/'profiles/corrected_profiles'),'--response',str(a.response),'--output',str(output/'calibration'),'--threads',str(a.threads),'--device',a.device]),
            ('permutations',[sys.executable,str(ROOT/'scripts/paper_permutations.py'),'--calibration',str(output/'calibration'),'--output',str(output/'permutations'),'--threads',str(a.threads),'--device',a.device])]
    results={}
    for name,command in stages:
        print('Running',name,flush=True)
        with (logs/(name+'.log')).open('w') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        results[name]=json.loads((output/name/'verification.json').read_text())
    paper.write_json(output/'verification.json',dict(status='passed',stages=results,frozen_model=True,fresh_training=False,fresh_calibration=True,fresh_permutations=True))
    print(output/'permutations/screening/reconstructed_candidates.tsv')

if __name__=='__main__':main()
