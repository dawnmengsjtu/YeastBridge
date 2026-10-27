#!/usr/bin/env python3
"""Assemble the six paper figures or redraw 42 data panels from complete templates."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.paper_release import digest,new_output,write_json
REL=Path('pd_product_figure_design/manuscript_review/yeastbridge_figures_20261010')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--redraw',action='store_true')
    p.add_argument('--rscript',default='Rscript')
    a=p.parse_args()
    spec=json.loads((ROOT/'data/paper/figure_assets.json').read_text())
    if digest(a.archive)!=spec['sha256']:p.error('Figure project archive checksum mismatch')
    import fitz
    output=new_output(a.output)
    with tarfile.open(a.archive,'r:gz') as archive:
        for item in archive.getmembers():
            path=output/item.name
            if not item.isfile() or not path.resolve().is_relative_to(output.resolve()):
                raise ValueError('Unsafe figure project member')
            path.parent.mkdir(parents=True,exist_ok=True)
            with archive.extractfile(item) as source,path.open('wb') as dest:shutil.copyfileobj(source,dest)
    project=output/'project';w=project/REL
    for row in json.loads((output/'manifest.json').read_text()):
        if digest(output/row['path'])!=row['sha256']:raise ValueError('Figure input checksum differs: '+row['path'])
    panels=json.loads((w/'panel_manifest.json').read_text())
    references=output/'native_reference';shutil.copytree(w/'build/expanded_native',references)
    logs=new_output(output/'logs')
    for sub in ['build/expanded_reflow','build/additions','build/main_native']:(w/sub).mkdir(parents=True,exist_ok=True)
    for sub in ['panels','pdf','preview']:(project/'delivery/yeastbridge/figures'/sub).mkdir(parents=True,exist_ok=True)
    def run(name,command,extra=None):
        print('Running',name,flush=True)
        with (logs/(name+'.log')).open('w') as log:
            subprocess.run(command,cwd=project,env=dict(os.environ,**(extra or {})),stdout=log,stderr=subprocess.STDOUT,check=True)
    if a.redraw:
        r_jobs=[('render_t5_update.R',None),('points_4_6.R','4b,4c,4d,4e,5b,5c,5d,6b,6c,6e,6f'),('candidate_points.R','candidate_docking_main,candidate_clusters,candidate_tail,candidate_bioactivity,candidate_barcodes,7b,7e,7f,7g,7h'),('go_points.R','go_b,go_c,go_d,go_e,go_f')]
        for script,wanted in r_jobs:
            run(script,[a.rscript,str(w/'source'/script)],dict(YB_USE_VERIFIED_SNAPSHOT='1',**({'YB_PANELS':wanted} if wanted else {})))
        run('evidence_wheel',[sys.executable,str(w/'source/evidence_main_wheel.py')])
        # The nine conceptual illustrations retain their approved native content.
        scenes={s['id']:s for s in json.loads((w/'workflow_scenes.json').read_text())}
        for id,row in panels.items():
            if not row['editable_workflow']:continue
            sid=row['source_id'];folder='build/additions' if sid.startswith('workflow_') else 'build/expanded_reflow'
            dest=w/folder/(sid+'.pdf');shutil.copy2(references/(id+'.pdf'),dest)
            if sid.startswith('workflow_'):write_json(dest.with_suffix('.scene.json'),scenes[id])
        run('integration',[sys.executable,str(w/'source/integrate_expansion.py')])
        # Integration reletters data panels; already-final schematic lettering is preserved.
        for id,row in panels.items():
            if row['editable_workflow']:shutil.copy2(references/(id+'.pdf'),w/'build/expanded_native'/(id+'.pdf'))
    run('assembly',[sys.executable,str(w/'source/deliver_expanded.py')],{'YB_PDF_ONLY':'1'})
    checks=[]
    for id,row in panels.items():
        with fitz.open(references/(id+'.pdf')) as old,fitz.open(w/'build/expanded_native'/(id+'.pdf')) as new:
            x=old[0].get_pixmap();y=new[0].get_pixmap()
            same=(x.width,x.height,x.samples)==(y.width,y.height,y.samples)
            checks.append(dict(panel=id,pixels_exact_at_72_dpi=same,recomputed_from_tables=bool(a.redraw and not row['editable_workflow'])))
    good=all(c['pixels_exact_at_72_dpi'] for c in checks)
    report=dict(status='passed' if good else 'rendering-difference',mode='redraw' if a.redraw else 'assemble',figures=6,panels=51,recomputed_data_panels=42 if a.redraw else 0,retained_schematic_panels=9 if a.redraw else 51,retained_molecular_render_assets=True,checks=checks,archive_sha256=spec['sha256'],script_sha256=digest(Path(__file__)))
    write_json(output/'verification.json',report)
    print(project/'delivery/yeastbridge/figures/pdf')
    return 0 if good else 2

if __name__=='__main__':raise SystemExit(main())
