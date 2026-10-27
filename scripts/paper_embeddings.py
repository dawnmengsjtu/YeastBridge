#!/usr/bin/env python3
"""Rebuild all frozen protein embeddings from the original sequence inputs."""
import argparse
import gzip
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import paper_release as paper


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',type=Path,required=True,help='Official ESM-2 650M checkpoint with adjacent contact-regression file')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--device',default='cuda')
    a=p.parse_args()
    paper.require('inference');paper.require('preprocessing')
    spec=json.loads((ROOT/'data/provenance/esm2/sequence-recovery.json').read_text())
    previous=json.loads((ROOT/'data/provenance/esm2/rebuild-run.json').read_text())
    if paper.digest(a.model)!=spec['esm2_checkpoint']['sha256']:p.error('ESM-2 checkpoint differs')
    regression=a.model.with_name(a.model.stem+'-contact-regression.pt')
    if not regression.is_file() or paper.digest(regression)!=previous['contact_regression_sha256']:p.error('ESM contact-regression file missing or changed')
    out=paper.new_output(a.output);inputs=paper.new_output(out/'inputs');logs=paper.new_output(out/'logs');checks=[]
    names={'human_810':'human_810.fasta.gz','universe_1175':'universe_1175.fasta.gz','yeast_650m':'yeast_6742.fasta.gz'}
    records={Path(row['file']).name:row for row in spec['files']}
    for name,filename in names.items():
        entry=records[filename];compressed=ROOT/entry['file']
        if paper.digest(compressed)!=entry['sha256']:raise ValueError('Frozen FASTA checksum differs')
        fasta=inputs/filename.removesuffix('.gz')
        with gzip.open(compressed,'rb') as source,fasta.open('wb') as dest:shutil.copyfileobj(source,dest)
        if paper.digest(fasta)!=entry['source_file_sha256']:raise ValueError('Decompressed FASTA checksum differs')
        command=[sys.executable,str(ROOT/'scripts/extract_esm2_candidates.py'),'--fasta',str(fasta),'--gene-master',str(ROOT/'raw/externalvalidation/mappings/gene_master.tsv'),'--model',str(a.model.resolve()),'--outdir',str(out/name),'--device',a.device]
        print('Rebuilding',name,flush=True)
        with (logs/(name+'.log')).open('w') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        for file in ['index.tsv','esm2_mean_fp32.npy']:
            reference=ROOT/'raw/tier1_esm2'/name/file
            sha=paper.digest(out/name/file);expected=paper.digest(reference)
            checks.append(dict(group=name,file=file,exact=sha==expected,reference_sha256=expected,rebuilt_sha256=sha))
    good=all(row['exact'] for row in checks)
    report=dict(status='passed' if good else 'different',proteins=8727,checks=checks,model_sha256=paper.digest(a.model),contact_regression_sha256=paper.digest(regression),packages={x:importlib.metadata.version(x) for x in ['torch','fair-esm','numpy','pandas','pyarrow']},script_sha256=paper.digest(Path(__file__)),extractor_sha256=paper.digest(ROOT/'scripts/extract_esm2_candidates.py'))
    paper.write_json(out/'verification.json',report);print(json.dumps(report,indent=2));return 0 if good else 2

if __name__=='__main__':raise SystemExit(main())
