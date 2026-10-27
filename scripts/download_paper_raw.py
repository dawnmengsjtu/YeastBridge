#!/usr/bin/env python3
"""Fetch public raw inputs and verify decompressed files against the paper snapshot."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import sys
import urllib.request
import zipfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.paper_release import digest,write_json


def fetch(url,path):
    pending=path.with_suffix(path.suffix+'.downloading')
    pending.parent.mkdir(parents=True,exist_ok=True)
    try:
        with urllib.request.urlopen(url,timeout=120) as response,pending.open('wb') as dest:
            shutil.copyfileobj(response,dest,1<<20)
        pending.replace(path)
    finally:pending.unlink(missing_ok=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',choices=['cel','counts'],required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-archives',type=int,help='Download only this many CEL archives for a partial access check')
    a=p.parse_args()
    if a.max_archives is not None and a.max_archives<1:p.error('--max-archives must be positive')
    spec=json.loads((ROOT/'data/paper/raw_sources.json').read_text())
    output=a.output.resolve();output.mkdir(parents=True,exist_ok=True)
    if a.dataset=='counts':
        entry=spec['counts'];target=output/entry['decompressed_name']
        if not target.exists():
            archive=output/(target.name+'.gz');fetch(entry['url'],archive)
            pending=target.with_suffix('.pending')
            try:
                with gzip.open(archive,'rb') as source,pending.open('wb') as dest:shutil.copyfileobj(source,dest)
                if digest(pending)!=entry['decompressed_sha256']:raise ValueError('Raw count checksum mismatch')
                pending.replace(target)
            finally:pending.unlink(missing_ok=True)
        if digest(target)!=entry['decompressed_sha256']:raise ValueError('Existing raw count file differs')
        report=dict(status='passed',dataset='GSE125162',files=1,sha256=digest(target))
    else:
        files={r['name']:r for r in spec['files']};done=set()
        cel=output/'cel_all';cel.mkdir(exist_ok=True)
        archives=output/'archives';archives.mkdir(exist_ok=True)
        selected=spec['archives'][:a.max_archives]
        for entry in selected:
            expected={n for n,r in files.items() if r['archive']==entry['name']}
            existing={n for n in expected if (cel/n).exists() and digest(cel/n)==files[n]['sha256']}
            if existing!=expected:
                archive=archives/entry['name']
                if not archive.exists():fetch(entry['url'],archive)
                with zipfile.ZipFile(archive) as z:
                    selected_members={}
                    for member in z.infolist():
                        name=Path(member.filename).name
                        if name in expected:
                            if name in selected_members:raise ValueError('Duplicate CEL basename')
                            selected_members[name]=member
                    if set(selected_members)!=expected:raise ValueError('Archive does not cover registered CEL files')
                    for name,member in selected_members.items():
                        if member.file_size!=files[name]['size']:raise ValueError('CEL size differs: '+name)
                        target=cel/name
                        if target.exists():
                            if digest(target)!=files[name]['sha256']:raise ValueError('Different existing CEL: '+name)
                            continue
                        pending=target.with_suffix('.pending')
                        try:
                            with z.open(member) as source,pending.open('wb') as dest:shutil.copyfileobj(source,dest)
                            if digest(pending)!=files[name]['sha256']:raise ValueError('CEL checksum differs: '+name)
                            pending.replace(target)
                        finally:pending.unlink(missing_ok=True)
            done.update(expected);print(entry['name'],len(done),'verified CEL files',flush=True)
        meta=output/'raw';meta.mkdir(exist_ok=True)
        src=ROOT/'raw/preprocessing/E-MTAB-2391.sdrf.txt';target=meta/src.name
        if target.exists() and digest(target)!=digest(src):raise ValueError('Different existing SDRF')
        if not target.exists():shutil.copyfile(src,target)
        report=dict(status='passed' if len(done)==len(files) else 'partial',dataset='E-MTAB-2391',verified_files=len(done),expected_files=len(files),archives=len(selected))
    write_json(output/'download_verification.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
