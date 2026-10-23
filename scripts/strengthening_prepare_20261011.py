from __future__ import annotations
import sys, json, hashlib, time, re, os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/public/home/mengxl/dzy/YeastBridge_0912')
UP = Path('/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901')
OUT = Path(os.environ.get('YB_ANALYSIS_DIR', ROOT/'results_a6000/strengthening_20261011'))
represented_only = os.environ.get('YB_REPRESENTED_ONLY') == '1'
OUT.mkdir(exist_ok=True)
if os.environ.get('YB_BINDINGDB_DOWNLOAD')=='1':
    import urllib.request
    name='BindingDB_All_202610_tsv.zip'
    base='https://www.bindingdb.org/rwd/bind/downloads/'
    path=OUT/name
    if not path.exists():
        partial=path.with_suffix('.zip.partial')
        with urllib.request.urlopen(base+name,timeout=60) as response,partial.open('wb') as handle:
            expected=int(response.headers.get('Content-Length','0'));total=0
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                handle.write(chunk);total+=len(chunk)
                if total%(32*1024*1024)==0:print('BindingDB bytes',total,flush=True)
        if expected and total!=expected:raise ValueError('incomplete BindingDB download')
        partial.rename(path)
    md5=hashlib.md5();sha256=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):md5.update(chunk);sha256.update(chunk)
    with urllib.request.urlopen(base+name.replace('.zip','.md5'),timeout=60) as response:source_md5=response.read().decode().strip().split()[0]
    if md5.hexdigest()!=source_md5:raise ValueError('BindingDB publisher checksum mismatch')
    report=dict(file=name,url=base+name,download_page='https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp',bytes=path.stat().st_size,sha256=sha256.hexdigest(),md5=md5.hexdigest(),publisher_md5=source_md5,retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    (OUT/'source_manifest.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
    sys.exit(0)
if os.environ.get('YB_BINDINGDB_PREPARE')=='1':
    import csv, io, zipfile
    candidates=pd.read_csv(os.environ['YB_CANDIDATE_CSV'],dtype=str,keep_default_na=False)
    candidate_keys=set(candidates.iloc[:,3]);candidate_cids=set(candidates.iloc[:,4])
    scope=pd.read_csv(os.environ['YB_BINDINGDB_COMPOUND_SCOPE'],sep='\t',dtype=str,keep_default_na=False)
    cohorts=scope.groupby('inchikey').cohort.agg(lambda s:'|'.join(sorted(set(s)))).to_dict()
    annotation=pd.read_csv(os.environ['YB_TARGET_ANNOTATIONS'],sep='\t',dtype=str,keep_default_na=False)
    annotation=annotation[annotation.screened.eq('True')]
    accession_targets={}
    for row in annotation.itertuples():
        for accession in row.uniprot_ids.split('|'):
            if accession:accession_targets.setdefault(accession,set()).add(row.target_id)
    candidate_rows=[];n_source=0;n_relevant=0;columns=None
    with zipfile.ZipFile(OUT/'BindingDB_All_202610_tsv.zip') as archive, (OUT/'library_measurements.tsv').open('w') as output:
        with archive.open('BindingDB_All.tsv') as member:
            reader=csv.reader(io.TextIOWrapper(member,encoding='utf-8'),delimiter='\t',quoting=csv.QUOTE_NONE)
            header=next(reader);lookup={name:i for i,name in enumerate(header)}
            key_i=lookup['Ligand InChI Key'];cid_i=lookup['PubChem CID']
            accession_i=[i for i,col in enumerate(header) if re.fullmatch(r'UniProt \((SwissProt|TrEMBL)\) Primary ID of Target Chain \d+',col)]
            kept=header[:52]
            columns=['source_row','query_cohort','candidate_exact_identity','candidate_cid_only','human_target_ids','single_human_target','target_annotation_flag']+kept
            writer=csv.DictWriter(output,fieldnames=columns,delimiter='\t');writer.writeheader()
            for n_source,values in enumerate(reader,1):
                if n_source%500000==0:print('BindingDB records',n_source,flush=True)
                if len(values)<=cid_i:raise ValueError('short BindingDB record at '+str(n_source))
                key=values[key_i].strip();cids=set(re.findall(r'\d+',values[cid_i]))
                exact=key in candidate_keys;cid_only=bool(cids & candidate_cids) and not exact
                if key not in cohorts and not exact and not cid_only:continue
                source=dict(zip(kept,values[:52]));accessions=set()
                for i in accession_i:
                    if i<len(values):accessions.update(re.split(r'[;,\s]+',values[i].strip()))
                targets=set().union(*(accession_targets.get(a,set()) for a in accessions))
                organism=source['Target Source Organism According to Curator or DataSource'].strip()
                single=source['Number of Protein Chains in Target (>1 implies a multichain complex)'].strip()=='1'
                flag='mutant_or_variant' if re.search(r'mutant|mutation|\[[^\]]*\]|\b[A-Z]\d+[A-Z]\b',source['Target Name'],re.I) else ''
                row=dict(source_row=n_source,query_cohort=cohorts.get(key,''),candidate_exact_identity=exact,candidate_cid_only=cid_only,human_target_ids='|'.join(sorted(targets)),single_human_target=single and len(targets)==1 and organism in ('','Homo sapiens'),target_annotation_flag=flag,**source)
                writer.writerow(row);n_relevant+=1
                if exact or cid_only:candidate_rows.append(row)
    measurements=pd.DataFrame(candidate_rows,columns=columns)
    measurements.to_csv(OUT/'candidate_compound_measurements.tsv',sep='\t',index=False)
    coverage=[];pair_records=[]
    for candidate in candidates.itertuples(index=False,name=None):
        candidate_id,target,key,cid=candidate[0],candidate[2],candidate[3],candidate[4]
        exact=measurements[measurements['Ligand InChI Key'].eq(key)]
        hit=exact.loc[exact.human_target_ids.map(lambda s:target in s.split('|')).astype(bool)]
        single=hit[hit.single_human_target & hit.target_annotation_flag.eq('')]
        cid_only=measurements[measurements['PubChem CID'].map(lambda s:cid in re.findall(r'\d+',s)) & measurements['Ligand InChI Key'].ne(key)]
        cid_target=cid_only.loc[cid_only.human_target_ids.map(lambda s:target in s.split('|')).astype(bool)]
        for row in hit.to_dict('records'):pair_records.append(dict(candidate_id=candidate_id,target_id=target,**row))
        coverage.append(dict(candidate_id=candidate_id,target_id=target,inchikey=key,cid=cid,n_exact_compound_records=len(exact),n_matching_target_records=len(hit),n_single_human_target_records=len(single),n_cid_only_target_records=len(cid_target),n_matching_source_documents=hit['Article DOI'].replace('',np.nan).nunique()))
    pd.DataFrame(pair_records,columns=['candidate_id','target_id']+columns).to_csv(OUT/'candidate_pair_measurements.tsv',sep='\t',index=False)
    coverage=pd.DataFrame(coverage);coverage.to_csv(OUT/'candidate_assay_coverage.tsv',sep='\t',index=False)
    summary=dict(n_source_records=n_source,n_library_measurements=n_relevant,n_scope_compounds=len(cohorts),n_nominees=len(candidates),n_candidate_compounds=len(candidate_keys),n_exact_candidate_compounds=measurements.loc[measurements.candidate_exact_identity,'Ligand InChI Key'].nunique(),n_candidate_pairs_with_matching_target=int(coverage.n_matching_target_records.gt(0).sum()),n_candidate_pairs_with_single_human_target=int(coverage.n_single_human_target_records.gt(0).sum()),model_scores_read=False,interpretation='Coverage audit only. BindingDB provenance must be deduplicated against original ChEMBL/PubChem source records before any independent validation claim.',inputs={name:dict(path=os.environ[name],sha256=hashlib.sha256(Path(os.environ[name]).read_bytes()).hexdigest()) for name in ['YB_CANDIDATE_CSV','YB_BINDINGDB_COMPOUND_SCOPE','YB_TARGET_ANNOTATIONS']})
    (OUT/'coverage_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
    sys.exit(0)
if os.environ.get('YB_GPCR_PANEL_PREPARE')=='1':
    from rdkit import Chem
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
    source=Path(os.environ['YB_GPCR_PANEL_SOURCE'])
    raw=pd.read_csv(source/'Screen_Summary_Pruned.csv',sep=';')
    scores=pd.read_csv(source/'Z-scores_summary_pruned.csv',sep=';')
    counter=pd.read_csv(source/'source_counter_screen.tsv',sep='\t',dtype=str)
    assert len(raw)==550 and len(scores)==539 and scores.catnr.is_unique
    names={'5HT4B':'HTR4','ADRA2B':'ADRA2B','OPRM1':'OPRM1','MTNR1A':'MTNR1A','ADRA2A':'ADRA2A','CHRM3':'CHRM3','5HT1A':'HTR1A'}
    matched=raw.set_index('catnr').loc[scores.catnr]
    errors={}
    for col in names:
        x=np.log10(matched[col].to_numpy(dtype=float));m=np.median(x)
        expected=.6745*(x-m)/np.median(abs(x-m))
        errors[col]=float(np.max(abs(expected-scores['Z'+col].to_numpy())))
    assert max(errors.values())<1e-9
    qc=[]
    for panel,a,b in [('A','J','K'),('B','F','G'),('C','A','B')]:
        frame=counter[[a,b]].rename(columns={a:'name',b:'rlu'}).dropna()
        frame['rlu']=pd.to_numeric(frame.rlu,errors='coerce');frame=frame.dropna()
        control=frame[frame.name.eq('Control')].rlu.mean()
        for name,g in frame.groupby('name'):
            if name=='Control':continue
            assert len(g)==3
            qc.append(dict(panel=panel,name=name,mean=g.rlu.mean(),sd=g.rlu.std(ddof=1),control_mean=control,excluded=g.rlu.mean()+g.rlu.std(ddof=1)<control))
    qc=pd.DataFrame(qc);assert qc.excluded.sum()==36
    norm=lambda s:re.sub('[^a-z0-9]','',s.lower())
    excluded=set()
    for i,row in qc[qc.excluded].iterrows():
        hit=scores[scores.Name.map(lambda s:norm(s).startswith(norm(row['name'])))]
        if row['name']=='Chlorpromazine':hit=hit[hit.Name.eq('Chlorpromazine_hydrochloride')]
        assert len(hit)==1,(row['name'],hit.Name.tolist())
        qc.loc[i,'catnr']=hit.catnr.iloc[0];excluded.add(hit.catnr.iloc[0])
    remover=SaltRemover.SaltRemover();uncharger=rdMolStandardize.Uncharger()
    def identity(smiles):
        mol=Chem.MolFromSmiles(smiles)
        if mol is None:raise ValueError('unparseable source structure')
        key=Chem.MolToInchiKey(mol);clean=remover.StripMol(mol,dontRemoveEverything=True)
        organic=[m for m in Chem.GetMolFrags(clean,asMols=True) if any(a.GetAtomicNum()==6 for a in m.GetAtoms())]
        parent=Chem.MolToInchiKey(uncharger.uncharge(rdMolStandardize.Cleanup(organic[0]))) if len(organic)==1 else ''
        return key,parent
    keys=[identity(s) for s in raw.Smiles];raw['inchikey']=[x[0] for x in keys];raw['parent_inchikey']=[x[1] for x in keys]
    duplicate_keys=set(raw.loc[raw.inchikey.duplicated(keep=False),'inchikey'])
    chemistry=raw.set_index('catnr').loc[scores.catnr].reset_index()
    chemistry['source_reporter_interference']=chemistry.catnr.isin(excluded)
    chemistry['source_duplicated_identity']=chemistry.inchikey.isin(duplicate_keys)
    rows=[]
    for i,record in enumerate(chemistry.to_dict('records')):
        for col,target in names.items():
            z=float(scores.iloc[i]['Z'+col]);eligible=not record['source_reporter_interference'] and not record['source_duplicated_identity'] and target!='HTR1A'
            rows.append(dict(catnr=record['catnr'],name=record['Name'],inchikey=record['inchikey'],parent_inchikey=record['parent_inchikey'],smiles=record['Smiles'],target_id=target,source_z=z,activation_hit=z>3,negative_response=z< -3,source_reporter_interference=record['source_reporter_interference'],source_duplicated_identity=record['source_duplicated_identity'],unreliable_source_receptor=target=='HTR1A',eligible=eligible))
    pd.DataFrame(rows).to_csv(OUT/'functional_panel_all_records.tsv',sep='\t',index=False)
    chemistry.to_csv(OUT/'functional_panel_compounds.tsv',sep='\t',index=False);qc.to_csv(OUT/'functional_panel_counter_screen.tsv',sep='\t',index=False)
    report=dict(source_catalog_rows=len(raw),source_unique_compounds=len(scores),source_qc_excluded_compounds=len(excluded),source_duplicate_identity_groups=len(duplicate_keys),eligible_compounds=int((~chemistry.source_reporter_interference & ~chemistry.source_duplicated_identity).sum()),eligible_receptors=6,max_formula_error=max(errors.values()),no_model_scores_read=True)
    (OUT/'functional_panel_preparation_summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
    sys.exit(0)
if os.environ.get('YB_PIOTROWSKI_PREPARE')=='1':
    import subprocess
    source=Path(os.environ['YB_PIOTROWSKI_SOURCE'])
    metadata=pd.read_excel(source/'dataset_4.xlsx',sheet_name='Suppl. dataset 4',header=5).fillna('')
    profiles=pd.read_excel(source/'dataset_5.xlsx',sheet_name='Suppl. dataset 5',header=2,index_col=0)
    ids=profiles.columns.str.split('_',n=1).str[0]
    if metadata.Condition_ID.duplicated().any() or ids.duplicated().any():raise ValueError('ambiguous external condition identity')
    metadata=metadata.set_index('Condition_ID').loc[ids].rename_axis('Condition_ID').reset_index()
    if not metadata.In_final_dataset.eq(True).all():raise ValueError('external source inclusion flag mismatch')
    orfs=profiles.index.str.rsplit('_',n=1).str[-1].to_numpy(dtype=str)
    if len(set(orfs))!=len(orfs):raise ValueError('duplicate external strain identity')
    complete=np.isfinite(profiles.to_numpy(dtype=float)).all(axis=1)
    missing=dict(zip(orfs[~complete],profiles.isna().sum(axis=1).to_numpy()[~complete].astype(int).tolist()))
    if profiles.shape!=(289,11983) or int(complete.sum())!=282 or sorted(missing.values())!=[3565]*7:raise ValueError('source missingness differs from frozen handling rule')
    # Structure preparation runs in the existing chemistry environment; human target labels and scores are not inputs.
    code='''import json,sys
from rdkit import Chem,rdBase
rows=[]
with rdBase.BlockLogs():
 for inchi,smiles in json.load(sys.stdin):
  a=Chem.MolFromInchi(inchi) if inchi else None
  b=Chem.MolFromSmiles(smiles) if smiles else None
  m=a if a is not None else b
  ka=Chem.MolToInchiKey(a) if a is not None else ''
  kb=Chem.MolToInchiKey(b) if b is not None else ''
  rows.append(dict(inchikey=Chem.MolToInchiKey(m) if m is not None else '',smiles=Chem.MolToSmiles(m) if m is not None else '',identity_basis='publisher_InChI' if a is not None else 'publisher_SMILES' if b is not None else 'unmapped',structure_disagreement=bool(ka and kb and ka!=kb),smiles_inchikey=kb))
print(json.dumps(rows))
'''
    identity=json.loads(subprocess.check_output([os.environ['YB_CHEMISTRY_PYTHON'],'-c',code],input=json.dumps(metadata[['InChI','SMILES']].values.tolist()).encode()))
    metadata=pd.concat([metadata,pd.DataFrame(identity)],axis=1)
    metadata['source_inchikey']=metadata.inchikey
    metadata.loc[metadata.structure_disagreement,'inchikey']=''
    values=profiles.to_numpy(dtype=float)[complete]
    if np.any(np.ptp(values,axis=0)==0):raise ValueError('constant external response: handling must be specified before scoring')
    vehicle=metadata.Name.str.strip().str.casefold().isin(['dmso','dmso control','negative control']).to_numpy()
    metadata.to_csv(OUT/'condition_metadata.tsv',sep='\t',index=False)
    metadata.loc[metadata.inchikey.ne(''),['inchikey','smiles']].drop_duplicates('inchikey').to_csv(OUT/'compounds.tsv',sep='\t',index=False)
    pd.DataFrame({'orf':orfs,'included':complete,'n_missing':profiles.isna().sum(axis=1).to_numpy()}).to_csv(OUT/'strain_coverage.tsv',sep='\t',index=False)
    np.savez_compressed(OUT/'strain_response.npz',strain_orfs=orfs[complete],compound_inchikeys=metadata.inchikey.to_numpy(dtype=str),doses=metadata.Concentration.to_numpy(dtype=str),array_names=ids.to_numpy(dtype=str),is_vehicle=vehicle,z_score=values.astype(np.float32),method=np.array(['Piotrowski2017 corrected publisher CG scores; 282 complete strains; no response re-normalization or imputation']))
    report=dict(n_source_strains=len(orfs),n_complete_strains=int(complete.sum()),omitted_strains=missing,n_conditions=len(metadata),n_unique_compounds=int(metadata.loc[metadata.inchikey.ne('')&~vehicle,'inchikey'].nunique()),n_unmapped=int(metadata.identity_basis.eq('unmapped').sum()),n_structure_disagreements=int(metadata.structure_disagreement.sum()),n_identity_excluded=int(metadata.inchikey.eq('').sum()),n_vehicle=int(vehicle.sum()),source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source/'dataset_2.xlsx',source/'dataset_4.xlsx',source/'dataset_5.xlsx',source/'analysis_plan.json']})
    report['n_unparseable_structures']=report['n_unmapped']
    report['n_unmapped']=int(metadata.source_inchikey.eq('').sum())
    report['n_parsed_without_full_identity']=report['n_unmapped']-report['n_unparseable_structures']
    (OUT/'source_preparation_summary.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True);sys.exit(0)
if os.environ.get('YB_HGNC_AUDIT')=='1':
    import urllib.request
    url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
    path=OUT/'hgnc_complete_set.txt'
    if not path.exists():
        with urllib.request.urlopen(url,timeout=120) as response:payload=response.read()
        if not payload.startswith(b'hgnc_id\tsymbol\t'):raise ValueError('unexpected HGNC response')
        path.write_bytes(payload)
    hgnc=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
    scope=pd.read_csv(UP/'product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t')
    rows=[]
    for target in scope.target_id:
        hits=hgnc[hgnc.symbol.eq(target)];basis='approved_symbol'
        if hits.empty:
            hits=hgnc[hgnc.prev_symbol.map(lambda s:target in s.split('|'))];basis='previous_symbol'
        if hits.empty:
            hits=hgnc[hgnc.alias_symbol.map(lambda s:target in s.split('|'))];basis='alias_symbol'
        row=dict(target_id=target,identity_basis=basis if len(hits)==1 else 'ambiguous' if len(hits)>1 else 'unmapped',n_matches=len(hits))
        for col in ['hgnc_id','symbol','name','locus_group','locus_type','status','entrez_id','uniprot_ids']:
            row[col]='|'.join(hits[col].astype(str))
        rows.append(row)
    result=scope.merge(pd.DataFrame(rows),on='target_id',validate='one_to_one')
    task_dir=Path(os.environ.get('YB_TASK_DIR',ROOT/'tasks/esm2_joint_tasks_dc_20260912'))
    screened={p.stem.replace('yeast_task_','') for p in task_dir.glob('yeast_task_*.tsv')}
    if not screened or not screened.issubset(set(result.target_id)):raise ValueError('task target scope mismatch')
    result['screened']=result.target_id.isin(screened)
    result.to_csv(OUT/'target_identity_annotations.tsv',sep='\t',index=False)
    current=result[result.screened]
    summary=dict(source_url=url,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),n_source_targets=len(result),n_screened_targets=len(current),not_screened=result.loc[~result.screened,'target_id'].tolist(),identity_counts=current.identity_basis.value_counts().to_dict(),locus_group_counts=current.locus_group.value_counts().to_dict(),scope_changed=False,interpretation='Current HGNC annotation only; no target removed and no outcome-based selection.')
    (OUT/'target_identity_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True);sys.exit(0)
if os.environ.get('YB_NIBR_DOWNLOAD')=='1':
    import urllib.request
    source=OUT/'source';source.mkdir(exist_ok=True)
    # Original 2014 release identifiers from the official Dryad landing page.
    files={'Table_S1.xls':67139,'HIP_scores.txt':67147,'HOP_scores.txt':67148}
    manifest=[]
    for name,file_id in files.items():
        path=source/name;url='https://datadryad.org/downloads/file_stream/'+str(file_id)
        if not path.exists():
            partial=path.with_suffix(path.suffix+'.partial')
            request=urllib.request.Request(url,headers={'User-Agent':'YeastBridge public-data replication'})
            with urllib.request.urlopen(request,timeout=120) as response,partial.open('wb') as handle:
                expected=int(response.headers.get('Content-Length','0'))
                total=0
                while True:
                    chunk=response.read(1024*1024)
                    if not chunk:break
                    handle.write(chunk);total+=len(chunk)
                    if total%(64*1024*1024)==0:print(name,total,'bytes',flush=True)
            if expected and total!=expected:raise ValueError('incomplete source download '+name)
            if total<10000:raise ValueError('unexpectedly small source download '+name)
            partial.rename(path)
        digest=hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda:handle.read(1024*1024),b''):digest.update(chunk)
        manifest.append(dict(file=name,url=url,bytes=path.stat().st_size,sha256=digest.hexdigest()))
        print('source verified',name,path.stat().st_size,flush=True)
    (OUT/'source_manifest.json').write_text(json.dumps({'doi':'10.5061/dryad.v5m8v','release':'2014-11-19','files':manifest},indent=2))
    sys.exit(0)
if os.environ.get('YB_PUBCHEM_AUDIT')=='1':
    import urllib.request,urllib.error,urllib.parse
    from concurrent.futures import ThreadPoolExecutor
    candidates=pd.read_csv(os.environ.get('YB_CANDIDATE_CSV',ROOT/'results.csv'))
    targets=sorted(candidates.iloc[:,2].unique());cids=sorted(set(candidates.iloc[:,4].astype(int)))
    base='https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
    cache=OUT/'pubchem';cache.mkdir(exist_ok=True)
    (cache/'analysis_plan.json').write_text(json.dumps({'scope':'all supplied frozen nominees; exact supplied CID and human gene identity','n_nominees':len(candidates),
        'source':'PubChem BioAssay compound assay summaries and human gene summaries',
        'interpretation':'all active, inactive and inconclusive records retained; primary screening hits and counterscreens require assay-level review before any target-engagement claim'},indent=2))
    def request(job):
        name,url=job;path=cache/(name+'.json')
        if path.exists():payload=path.read_bytes()
        else:
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'YeastBridge scientific audit'}),timeout=60) as r:payload=r.read()
                    json.loads(payload);path.write_bytes(payload);break
                except urllib.error.HTTPError as e:
                    if e.code==404:payload=e.read();path.write_bytes(payload);break
                    if attempt==2:raise
                    time.sleep(2*(attempt+1))
                except Exception:
                    if attempt==2:raise
                    time.sleep(2*(attempt+1))
        return name,json.loads(payload),dict(url=url,path=str(path),sha256=hashlib.sha256(payload).hexdigest())
    jobs=[('gene_'+t,base+'gene/genesymbol/'+urllib.parse.quote(t)+'/9606/summary/JSON') for t in targets]
    jobs += [('cid_'+str(c),base+'compound/cid/'+str(c)+'/assaysummary/JSON') for c in cids]
    records={};manifest=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for i,(name,data,entry) in enumerate(pool.map(request,jobs),1):
            records[name]=data;manifest.append(entry)
            if i%10==0:print('PubChem snapshots',i,len(jobs),flush=True)
    genes=[];tables=[]
    for target in targets:
        hits=[g for g in records['gene_'+target].get('GeneSummaries',{}).get('GeneSummary',[]) if g.get('TaxonomyID')==9606 and (g.get('Symbol')==target or target in g.get('Synonym',[]))]
        if len(hits)>1:raise ValueError('ambiguous human gene '+target)
        genes.append(dict(target_id=target,gene_id=hits[0]['GeneID'] if hits else '',current_symbol=hits[0]['Symbol'] if hits else '',
            gene_name=hits[0].get('Name','') if hits else '',pseudogene_annotation='pseudogene' in hits[0].get('Name','').lower() if hits else None))
    for cid in cids:
        tab=records['cid_'+str(cid)].get('Table')
        if tab:
            frame=pd.DataFrame([r['Cell'] for r in tab['Row']],columns=tab['Columns']['Column']);frame['query_cid']=cid;tables.append(frame)
    all_rows=pd.concat(tables,ignore_index=True);all_rows.to_csv(cache/'all_assay_records.tsv',sep='\t',index=False)
    gene_map={r['target_id']:str(r['gene_id']) for r in genes};matches=[];counts=[]
    for r in candidates.itertuples(index=False,name=None):
        target,cid=r[2],int(r[4]);g=all_rows[all_rows.query_cid.eq(cid)]
        gene=gene_map[target];hit=g.loc[g['Target GeneID'].map(lambda x:gene in re.findall(r'\d+',str(x))).astype(bool)] if gene else g.iloc[:0]
        for row in hit.to_dict('records'):matches.append(dict(candidate_id=r[0],target_id=target,**row))
        counts.append(dict(candidate_id=r[0],target_id=target,cid=cid,human_gene_id=gene,n_assay_rows=len(hit),n_active_assays=hit[hit['Activity Outcome'].eq('Active')].AID.nunique(),n_inactive_assays=hit[hit['Activity Outcome'].eq('Inactive')].AID.nunique()))
    pd.DataFrame(genes).to_csv(cache/'human_gene_mapping.tsv',sep='\t',index=False)
    annotation=pd.DataFrame([dict(candidate_id=r[0],target_id=r[2]) for r in candidates.itertuples(index=False,name=None)]).merge(pd.DataFrame(genes),on='target_id',validate='many_to_one')
    annotation['source_url']=annotation.gene_id.map(lambda x:'https://www.ncbi.nlm.nih.gov/gene/'+str(x) if x else '')
    annotation.to_csv(cache/'candidate_gene_annotations.tsv',sep='\t',index=False)
    pd.DataFrame(matches,columns=['candidate_id','target_id']+all_rows.columns.tolist()).to_csv(cache/'candidate_assay_records.tsv',sep='\t',index=False)
    pd.DataFrame(counts).to_csv(cache/'candidate_assay_coverage.tsv',sep='\t',index=False)
    (cache/'snapshot_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(pd.DataFrame(counts).query('n_assay_rows>0').to_string(index=False),flush=True);sys.exit(0)
if os.environ.get('YB_REFERENCE_LIGAND_CONTEXT')=='1':
    from rdkit import Chem,DataStructs,rdBase
    from rdkit.Chem import SaltRemover,MACCSkeys,rdFingerprintGenerator
    from rdkit.Chem.MolStandardize import rdMolStandardize
    remover=SaltRemover.SaltRemover();uncharger=rdMolStandardize.Uncharger()
    def parent_mol(smiles):
        mol=Chem.MolFromSmiles(smiles or '')
        if mol is None or mol.GetNumAtoms()==0:return None
        mol=remover.StripMol(mol,dontRemoveEverything=True)
        parts=Chem.GetMolFrags(mol,asMols=True,sanitizeFrags=True)
        if len(parts)!=1:return None
        return uncharger.uncharge(rdMolStandardize.Cleanup(parts[0]))
    assays={r['assay_chembl_id']:r for r in json.loads((OUT/'reference_assays.json').read_text())}
    mapping=pd.read_csv(OUT/'reference_target_mapping.tsv',sep='\t').set_index('target_chembl_id').target_id.to_dict()
    records=[];mol_cache={}
    for r in json.loads((OUT/'reference_activities.json').read_text()):
        assay=assays[r['assay_chembl_id']];reasons=[]
        if assay.get('confidence_score')!=9:reasons.append('not_direct_confidence_9')
        if assay.get('assay_organism') not in [None,'','Homo sapiens']:reasons.append('nonhuman_assay')
        if r.get('assay_variant_mutation') or assay.get('variant_sequence'):reasons.append('variant')
        if r.get('standard_flag')!=1 or r.get('potential_duplicate'):reasons.append('nonstandard_or_duplicate')
        if r.get('data_validity_comment') not in [None,'','Manually validated']:reasons.append('validity_flag')
        value=float(r['standard_value'])
        if not 0<value<=1000 or r.get('standard_relation') not in ['=','<','<=']:reasons.append('not_potent_value_or_upper_bound')
        if 'inactive' in (r.get('activity_comment') or '').lower():reasons.append('inactive_comment')
        mid=r['molecule_chembl_id']
        if mid not in mol_cache:mol_cache[mid]=parent_mol(r.get('canonical_smiles'))
        mol=mol_cache[mid]
        if mol is None:reasons.append('ambiguous_or_missing_structure')
        records.append(dict(target_id=mapping[r['target_chembl_id']],molecule_chembl_id=mid,assay_type=r['assay_type'],activity_id=r['activity_id'],assay_chembl_id=r['assay_chembl_id'],document_chembl_id=r['document_chembl_id'],standard_type=r['standard_type'],standard_value_nM=value,standard_relation=r['standard_relation'],eligible=not reasons,exclusion_reasons='|'.join(reasons),parent_inchikey=Chem.MolToInchiKey(mol) if mol is not None else ''))
    frame=pd.DataFrame(records);frame.to_csv(OUT/'reference_ligand_curation.tsv',sep='\t',index=False)
    accepted=frame[frame.eligible];references=accepted.drop_duplicates(['target_id','assay_type','parent_inchikey'])
    chemistry=pd.read_csv(OUT.parent/'parent_identity_benchmark/screen_compound_parent_identity.tsv',sep='\t',keep_default_na=False)
    library=[parent_mol(s) for s in chemistry.smiles];valid=np.array([m is not None for m in library]);keys=chemistry.inchikey.to_numpy();cp={k:i for i,k in enumerate(keys)}
    nominated=pd.read_csv(OUT.parent/'candidate_score_decomposition.tsv',sep='\t')[['direction','target_id','inchikey']]
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    fps={'Morgan_r2_2048':generator.GetFingerprint,'MACCS_166':MACCSkeys.GenMACCSKeys}
    all_scores={};candidate_rows=[]
    for fingerprint,fp in fps.items():
        library_fp=[fp(m) if m is not None else None for m in library]
        for (target,assay_type),g in references.groupby(['target_id','assay_type']):
            ref=g.to_dict('records');ref_fps=[fp(mol_cache[r['molecule_chembl_id']]) for r in ref]
            maxima=np.full(len(keys),np.nan);best=np.full(len(keys),-1)
            for i in np.flatnonzero(valid):
                similarities=np.array(DataStructs.BulkTanimotoSimilarity(library_fp[i],ref_fps));best[i]=int(np.argmax(similarities));maxima[i]=similarities[best[i]]
            all_scores[fingerprint+'|'+assay_type+'|'+target]=maxima
            for r in nominated[nominated.target_id.eq(target)].itertuples():
                i=cp[r.inchikey];sim=maxima[i];nearest=ref[best[i]] if valid[i] else {}
                candidate_rows.append(dict(direction=r.direction,target_id=target,inchikey=r.inchikey,assay_type=assay_type,fingerprint=fingerprint,n_reference_ligands=len(ref),max_tanimoto=sim,library_percentile=float((np.sum(maxima[valid]<sim)+.5*np.sum(maxima[valid]==sim))/valid.sum()) if valid[i] else np.nan,nearest_molecule=nearest.get('molecule_chembl_id',''),nearest_parent_key=nearest.get('parent_inchikey',''),nearest_activity_id=nearest.get('activity_id',''),nearest_assay_id=nearest.get('assay_chembl_id',''),nearest_activity_nM=nearest.get('standard_value_nM',''),valid_structure=bool(valid[i])))
    np.savez_compressed(OUT/'reference_ligand_similarity.npz',inchikeys=keys,**all_scores)
    result=pd.DataFrame(candidate_rows)
    membership=pd.read_csv(OUT.parent/'reproduced_profile_family_membership.tsv',sep='\t').set_index('target_id').profile_family
    for i,r in result.iterrows():
        siblings=membership[membership.eq(membership.loc[r.target_id])].index
        alternatives=[all_scores[r.fingerprint+'|'+r.assay_type+'|'+t][cp[r.inchikey]] for t in siblings if t!=r.target_id and r.fingerprint+'|'+r.assay_type+'|'+t in all_scores]
        result.loc[i,'n_reference_covered_sibling_targets']=len(alternatives)
        result.loc[i,'margin_over_best_covered_sibling']=r.max_tanimoto-max(alternatives) if alternatives else np.nan
    result.to_csv(OUT/'candidate_reference_ligand_context.tsv',sep='\t',index=False)
    nominated.merge(result,on=['direction','target_id','inchikey'],how='left').to_csv(OUT/'all_nominees_ligand_coverage.tsv',sep='\t',index=False)
    summary=dict(rdkit=rdBase.rdkitVersion,n_downloaded_records=len(frame),n_eligible_records=len(accepted),n_reference_target_assay_pairs=references[['target_id','assay_type']].drop_duplicates().shape[0],n_distinct_reference_ligands=references.parent_inchikey.nunique(),n_nominated_pairs_with_reference_ligands=result[['target_id','inchikey']].drop_duplicates().shape[0],n_valid_library_structures=int(valid.sum()),interpretation='Descriptive ligand-chemistry support; library percentile is not a significance P value; activity of the nominated compound is not inferred from analog potency. All nominees remain in the coverage table.')
    (OUT/'reference_ligand_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True);print(result.sort_values('max_tanimoto',ascending=False).head(20).to_string(index=False),flush=True);sys.exit(0)
if os.environ.get('YB_CHEMBL_CURATE') == '1':
    # Freeze experimental evidence without loading any YeastBridge scores.
    plan={'source':'ChEMBL 37 public API snapshot', 'identity':'Original or conservative parent full InChIKey; ChEMBL parent hierarchy, never stereo-parent hierarchy',
          'target':'Human SINGLE PROTEIN; exact gene symbol in frozen 1175-target scope; direct assay confidence 9; no annotated variant',
          'activity':'Binding and functional assays retained separately; Ki/Kd/IC50/EC50 in nM; valid positive numeric values; relation respected; potential duplicate records excluded',
          'classification':'Strong active: <=1000 nM; measured weak/inactive: >=10000 nM; any active plus weak evidence is conflicting; intermediate/ambiguous records remain unclassified',
          'score_selection':'No screen score is read during evidence curation; original screen and benchmark remain unchanged'}
    (OUT/'curation_plan.json').write_text(json.dumps(plan,indent=2))
    scope=set(pd.read_csv(UP/'product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t').target_id)
    if os.environ.get('YB_EXTERNAL_PROTOCOL'):
        scope={p.stem.replace('yeast_task_','') for p in Path(os.environ['YB_TASK_DIR']).glob('yeast_task_*.tsv')}
    target_map={}
    for r in json.loads((OUT/'targets.json').read_text()):
        symbols={s['component_synonym'] for c in r.get('target_components',[]) for s in c.get('target_component_synonyms',[]) if s.get('syn_type')=='GENE_SYMBOL'}
        hit=symbols & scope
        if r.get('organism')=='Homo sapiens' and r.get('target_type')=='SINGLE PROTEIN' and len(hit)==1:
            target_map[r['target_chembl_id']]=next(iter(hit))
    assays={r['assay_chembl_id']:r for r in json.loads((OUT/'assays.json').read_text())}
    molecules=json.loads((OUT/'molecules.json').read_text())
    molecule_by_key={r['molecule_structures']['standard_inchi_key']:r for r in molecules}
    chemistry=pd.read_csv(OUT.parent/'parent_identity_benchmark/screen_compound_parent_identity.tsv',sep='\t',keep_default_na=False)
    identity=[]
    for r in chemistry.itertuples():
        for key,basis in [(r.inchikey,'original_exact'),(r.parent_inchikey,'standardized_parent')]:
            if key not in molecule_by_key:continue
            m=molecule_by_key[key];hierarchy=m.get('molecule_hierarchy') or {}
            identity.append({'inchikey':r.inchikey,'molecule_chembl_id':m['molecule_chembl_id'],
                             'parent_chembl_id':hierarchy.get('parent_chembl_id') or m['molecule_chembl_id'],
                             'identity_basis':basis,'compound_name':m.get('pref_name') or ''})
    identity=pd.DataFrame(identity).drop_duplicates();identity.to_csv(OUT/'chemical_identity.tsv',sep='\t',index=False)
    parent_map=identity.groupby('parent_chembl_id').inchikey.agg(lambda x:sorted(set(x))).to_dict()
    rows=[];seen=set()
    for r in json.loads((OUT/'activities.json').read_text()):
        if r['activity_id'] in seen:continue
        seen.add(r['activity_id'])
        if r['target_chembl_id'] not in target_map:continue
        keys=parent_map.get(r.get('parent_molecule_chembl_id') or r['molecule_chembl_id'],[])
        if not keys:continue
        assay=assays.get(r['assay_chembl_id'],{});reasons=[]
        if assay.get('confidence_score')!=9:reasons.append('not_direct_confidence_9')
        if assay.get('assay_organism') not in [None,'','Homo sapiens']:reasons.append('nonhuman_assay')
        if r.get('assay_variant_mutation') or assay.get('variant_sequence'):reasons.append('variant')
        if r.get('standard_flag')!=1:reasons.append('not_standardized')
        if r.get('potential_duplicate'):reasons.append('potential_duplicate')
        if r.get('data_validity_comment') not in [None,'','Manually validated']:reasons.append('validity_flag')
        try:value=float(r['standard_value'])
        except (ValueError,TypeError):value=np.nan
        relation=r.get('standard_relation')
        if r.get('standard_units')!='nM' or not np.isfinite(value) or value<=0:reasons.append('not_numeric_nM')
        if relation not in ['=','<','<=','>','>=']:reasons.append('ambiguous_relation')
        comment=(r.get('activity_comment') or '').lower()
        active=value<=1000 and relation in ['=','<','<=']
        weak=value>=10000 and relation in ['=','>','>=']
        if 'inactive' in comment and active:reasons.append('comment_numeric_conflict')
        for key in keys:
            rows.append(dict(inchikey=key,target_id=target_map[r['target_chembl_id']],activity_id=r['activity_id'],
                             assay_chembl_id=r['assay_chembl_id'],target_chembl_id=r['target_chembl_id'],molecule_chembl_id=r['molecule_chembl_id'],
                             document_chembl_id=r['document_chembl_id'],assay_type=r['assay_type'],standard_type=r['standard_type'],
                             standard_relation=relation,standard_value_nM=value,activity_comment=r.get('activity_comment'),
                             eligible=not reasons,exclusion_reasons='|'.join(reasons),active=active,weak=weak,
                             assay_description=r['assay_description'],document_year=r['document_year']))
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'curation_all_records.tsv',sep='\t',index=False)
    accepted=frame[frame.eligible].copy();accepted.to_csv(OUT/'curated_activity_records.tsv',sep='\t',index=False)
    summaries=[]
    for (key,target,typ),g in accepted.groupby(['inchikey','target_id','assay_type']):
        active=bool(g.active.any());weak=bool(g.weak.any())
        label='conflicting' if active and weak else 'active' if active else 'weak' if weak else 'intermediate'
        summaries.append({'inchikey':key,'target_id':target,'assay_type':typ,'evidence_class':label,
                          'n_records':len(g),'n_documents':g.document_chembl_id.nunique(),
                          'activity_ids':'|'.join(g.activity_id.astype(str)),'document_ids':'|'.join(sorted(set(g.document_chembl_id)))})
    pairs=pd.DataFrame(summaries);pairs.to_csv(OUT/'experimental_pair_evidence.tsv',sep='\t',index=False)
    if not os.environ.get('YB_EXTERNAL_PROTOCOL') or (OUT.parent/'candidate_score_decomposition.tsv').exists():
        candidates=pd.read_csv(OUT.parent/'candidate_score_decomposition.tsv',sep='\t')[['direction','target_id','inchikey']]
        candidates.merge(pairs,on=['inchikey','target_id'],how='left').to_csv(OUT/'candidate_experimental_evidence.tsv',sep='\t',index=False)
    summary={'n_downloaded_activity_records':len(seen),'n_in_scope_records':len(frame),'n_eligible_records':len(accepted),
             'n_eligible_compounds':accepted.inchikey.nunique(),'n_eligible_targets':accepted.target_id.nunique(),
             'pair_classes':pairs.groupby(['assay_type','evidence_class']).size().to_dict()}
    summary['pair_classes']={str(k):int(v) for k,v in summary['pair_classes'].items()}
    (OUT/'curation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
    sys.exit(0)
if os.environ.get('YB_CHEMBL_SNAPSHOT'):
    import urllib.request, urllib.parse
    from concurrent.futures import ThreadPoolExecutor
    stage=os.environ['YB_CHEMBL_SNAPSHOT']
    base='https://www.ebi.ac.uk/chembl/api/data/'
    cache=OUT/stage;cache.mkdir(exist_ok=True)
    def fetch(job):
        name,url=job;path=cache/(name+'.json')
        if path.exists():payload=path.read_bytes()
        else:
            for attempt in range(3):
                try:
                    request=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'YeastBridge scientific validation'})
                    with urllib.request.urlopen(request,timeout=45) as response:payload=response.read()
                    json.loads(payload);path.write_bytes(payload);break
                except Exception:
                    if attempt==2:raise
                    time.sleep(2*(attempt+1))
        return name,json.loads(payload),{'url':url,'path':str(path),'sha256':hashlib.sha256(payload).hexdigest()}
    if stage in ['human_targets','reference_activities']:
        records=[];manifest=[]
        if stage=='human_targets':
            urls=[base+'target.json?'+urllib.parse.urlencode({'organism':'Homo sapiens','target_type':'SINGLE PROTEIN','limit':1000})];field='targets'
        else:
            nominated=set(pd.read_csv(OUT.parent/'candidate_score_decomposition.tsv',sep='\t').target_id)
            membership=pd.read_csv(OUT.parent/'reproduced_profile_family_membership.tsv',sep='\t')
            families=set(membership[membership.target_id.isin(nominated)].profile_family)
            scope=set(membership[membership.profile_family.isin(families)].target_id)
            mapping=[]
            for r in json.loads((OUT/'human_targets.json').read_text()):
                symbols={v['component_synonym'] for c in r.get('target_components',[]) for v in c.get('target_component_synonyms',[]) if v.get('syn_type')=='GENE_SYMBOL'}
                hit=symbols & scope
                if len(hit)==1:mapping.append({'target_id':next(iter(hit)),'target_chembl_id':r['target_chembl_id']})
            table=pd.DataFrame(mapping).drop_duplicates();table.to_csv(OUT/'reference_target_mapping.tsv',sep='\t',index=False)
            ids=sorted(set(table.target_chembl_id));field='activities'
            plan={'scope':'all 60 frozen nominated targets and all members of their frozen profile clusters','reference':'human single-protein targets; binding and functional Ki/Kd/IC50/EC50 at <=1000 nM, assay confidence and record validity checked after retrieval','purpose':'independent ligand-chemistry context, not experimental validation of nominees','nomination':'original 80 pairs unchanged','source':'10.1038/nature08506 provides a ligand-comparison precedent'}
            (OUT/'reference_ligand_plan.json').write_text(json.dumps(plan,indent=2))
            urls=[base+'activity.json?'+urllib.parse.urlencode({'target_chembl_id__in':','.join(ids[i:i+10]),'target_organism':'Homo sapiens','assay_type__in':'B,F','standard_type__in':'Ki,Kd,IC50,EC50','standard_units':'nM','standard_value__lte':1000,'limit':1000}) for i in range(0,len(ids),10)]
        for batch,url in enumerate(urls):
            page=0
            while url:
                _,data,record=fetch((f'batch_{batch:03d}_page_{page:03d}',url));page+=1
                records.extend(data[field]);manifest.append(record)
                nxt=data['page_meta']['next'];url=urllib.parse.urljoin(base,nxt) if nxt else None
                print(stage,'batch',batch,'page',page,'records',len(records),flush=True)
        (OUT/(stage+'.json')).write_text(json.dumps(records))
    elif stage=='molecules':
        chemistry=pd.read_csv(OUT.parent/'parent_identity_benchmark/screen_compound_parent_identity.tsv',sep='\t',keep_default_na=False)
        keys=sorted({v for c in ['inchikey','parent_inchikey'] for v in chemistry[c] if v})
        jobs=[(f'batch_{i//40:03d}',base+'molecule.json?'+urllib.parse.urlencode({'molecule_structures__standard_inchi_key__in':','.join(keys[i:i+40]),'limit':1000})) for i in range(0,len(keys),40)]
        records=[];manifest=[]
        with ThreadPoolExecutor(max_workers=3) as pool:
            for i,(name,data,record) in enumerate(pool.map(fetch,jobs),1):
                if data['page_meta']['next']:raise ValueError('unexpected molecule pagination')
                records.extend(data['molecules']);manifest.append(record)
                if i%10==0:print('ChEMBL molecule batches',i,len(jobs),flush=True)
        by_id={r['molecule_chembl_id']:r for r in records}
        (OUT/'molecules.json').write_text(json.dumps(list(by_id.values())))
    elif stage=='activities':
        records=json.loads((OUT/'molecules.json').read_text())
        ids=sorted({r['molecule_hierarchy']['parent_chembl_id'] for r in records if r.get('molecule_hierarchy') and r['molecule_hierarchy'].get('parent_chembl_id')})
        orphan=sorted({r['molecule_chembl_id'] for r in records if not r.get('molecule_hierarchy') or not r['molecule_hierarchy'].get('parent_chembl_id')})
        batches=[{field:','.join(values[i:i+20])} for field,values in [('parent_molecule_chembl_id__in',ids),('molecule_chembl_id__in',orphan)] for i in range(0,len(values),20)]
        manifest=[];activities=[]
        def activity_batch(item):
            batch,identity=item;batch_records=[];batch_manifest=[]
            i=batch*20
            url=base+'activity.json?'+urllib.parse.urlencode(dict(identity,target_organism='Homo sapiens',assay_type__in='B,F',standard_type__in='Ki,Kd,IC50,EC50',limit=1000))
            page=0
            while url:
                name,data,record=fetch((f'batch_{i//20:03d}_page_{page:03d}',url));page+=1
                batch_records.extend(data['activities']);batch_manifest.append(record)
                nxt=data['page_meta']['next'];url=urllib.parse.urljoin(base,nxt) if nxt else None
            return batch_records,batch_manifest
        with ThreadPoolExecutor(max_workers=3) as pool:
            for batch,(batch_records,batch_manifest) in enumerate(pool.map(activity_batch,enumerate(batches))):
                activities.extend(batch_records);manifest.extend(batch_manifest)
                print('ChEMBL activities batch',batch,'total records',len(activities),flush=True)
        (OUT/'activities.json').write_text(json.dumps(activities))
    elif stage in ['targets','assays','reference_assays']:
        activities=json.loads((OUT/('reference_activities.json' if stage=='reference_assays' else 'activities.json')).read_text())
        field='target_chembl_id' if stage=='targets' else 'assay_chembl_id'
        endpoint='target' if stage=='targets' else 'assay'
        if stage=='assays':
            scope=set(pd.read_csv(UP/'product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t').target_id)
            target_records=json.loads((OUT/'targets.json').read_text())
            accepted=set()
            for r in target_records:
                symbols={s['component_synonym'] for c in r.get('target_components',[]) for s in c.get('target_component_synonyms',[]) if s.get('syn_type')=='GENE_SYMBOL'}
                if r.get('organism')=='Homo sapiens' and r.get('target_type')=='SINGLE PROTEIN' and len(symbols & scope)==1:
                    accepted.add(r['target_chembl_id'])
            if not accepted:raise ValueError('no ChEMBL targets mapped: verify target metadata fields')
            activities=[r for r in activities if r['target_chembl_id'] in accepted]
        ids=sorted({r[field] for r in activities});records=[];manifest=[]
        jobs=[(f'batch_{i//40:03d}',base+endpoint+'.json?'+urllib.parse.urlencode({field+'__in':','.join(ids[i:i+40]),'limit':1000})) for i in range(0,len(ids),40)]
        with ThreadPoolExecutor(max_workers=3) as pool:
            for name,data,record in pool.map(fetch,jobs):
                if data['page_meta']['next']:raise ValueError('unexpected metadata pagination')
                records.extend(data['assays' if stage=='reference_assays' else stage]);manifest.append(record)
        (OUT/(stage+'.json')).write_text(json.dumps(records))
    else:raise ValueError('unsupported ChEMBL snapshot stage')
    (OUT/(stage+'_manifest.json')).write_text(json.dumps(manifest,indent=2))
    print('ChEMBL snapshot complete',stage,len(activities) if stage=='activities' else len(records),flush=True)
    sys.exit(0)
if os.environ.get('YB_BENCHMARK_PARENT') == '1':
    from rdkit import Chem, rdBase
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
    # Chemical identity is frozen before any target-ranking scores are read.
    md=UP/'drug_ko_benchmark_v1/source_metadata'
    sample=pd.read_csv(md/'repo-sample-annotation-20250818.txt',sep='\t',comment='!',dtype=str).fillna('')
    drug=pd.read_csv(md/'repo-drug-annotation-20250818.txt',sep='\t',comment='!',dtype=str).fillna('')
    chemistry=pd.read_csv(os.environ.get('YB_COMPOUND_TABLE','/public/home/mengxl/dzy/yeastbridge_vs/data/external/paired_compounds/hiphop_lee2014_e_mtab_2391/processed/compounds.tsv.gz'),sep='\t',dtype=str).fillna('')
    chemistry=chemistry[['inchikey','smiles']].drop_duplicates()
    chemistry=chemistry[chemistry.inchikey.str.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]')].copy()
    if not os.environ.get('YB_EXTERNAL_PROTOCOL') and len(chemistry)!=3250:raise ValueError('screen compound identity scope differs from 3250')
    if chemistry.inchikey.duplicated().any():raise ValueError('ambiguous source compound structure')
    scope=pd.read_csv(UP/'product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t')
    targets=set(scope.target_id)
    if os.environ.get('YB_EXTERNAL_PROTOCOL'):
        targets={p.stem.replace('yeast_task_','') for p in Path(os.environ['YB_TASK_DIR']).glob('yeast_task_*.tsv')}
        if len(targets)!=1175:raise ValueError('external evaluation changed human target universe')
    remover=SaltRemover.SaltRemover();uncharger=rdMolStandardize.Uncharger()
    def parent(s):
        mol=Chem.MolFromSmiles(s)
        if mol is None:return '', 'invalid'
        mol=remover.StripMol(mol,dontRemoveEverything=True)
        fragments=Chem.GetMolFrags(mol,asMols=True)
        organic=[m for m in fragments if any(a.GetAtomicNum()==6 for a in m.GetAtoms())]
        if len(organic)!=1:return '', 'ambiguous_or_nonorganic'
        mol=uncharger.uncharge(rdMolStandardize.Cleanup(organic[0]))
        return Chem.MolToInchiKey(mol), 'ok'
    for frame,col in [(chemistry,'smiles'),(sample,'smiles')]:
        normalized=[parent(x) for x in frame[col]]
        frame['parent_inchikey']=[x[0] for x in normalized]
        frame['identity_status']=[x[1] for x in normalized]
    chemistry.to_csv(OUT/'screen_compound_parent_identity.tsv',sep='\t',index=False)
    sample.to_csv(OUT/'hub_parent_identity.tsv',sep='\t',index=False)
    matched=chemistry[(chemistry.identity_status=='ok') & chemistry.parent_inchikey.ne('')].merge(sample[(sample.identity_status=='ok') & sample.parent_inchikey.ne('')],on='parent_inchikey',suffixes=('_screen','_hub'))
    matched['identity_basis']='standardized_parent'
    exact=chemistry.merge(sample,left_on='inchikey',right_on='InChIKey',suffixes=('_screen','_hub'))
    exact['identity_basis']='exact_original_InChIKey'
    matched=pd.concat([exact,matched],ignore_index=True).merge(drug,on='pert_iname',how='left').fillna('')
    matched.to_csv(OUT/'benchmark_annotation_overlap.tsv',sep='\t',index=False)
    rows=[]
    for key,frame in matched.groupby('inchikey'):
        gold=sorted({g.strip() for v in frame.target for g in re.split(r'[|;,]',v) if g.strip() in targets})
        allgold=sorted({g.strip() for v in frame.target for g in re.split(r'[|;,]',v) if g.strip()})
        for g in gold:rows.append({'inchikey':key,'target_id':g,'compound_names':'|'.join(sorted(set(frame.pert_iname))),'moa':'|'.join(sorted(set(frame.moa))),'all_annotated_targets':'|'.join(allgold),'in_scope_gold_count':len(gold)})
    gold=pd.DataFrame(rows);gold.to_csv(OUT/'benchmark_frozen_pairs.tsv',sep='\t',index=False)
    plan={'question':'Does conservative parent-compound matching recover omitted external annotations?',
          'identity':'Union of original exact InChIKey and conservative parent standard InChIKey; defined-salt removal; reject ambiguous multiple organic parents; preserve stereochemistry and isotopes; standard InChI handles supported tautomer equivalence',
          'rdkit':rdBase.rdkitVersion,'n_compounds':int(gold.inchikey.nunique()),'n_pairs':len(gold),
          'rank_endpoint':'unchanged bidirectional maximum; rank all 1175 targets; all four prespecified metrics and directions retained',
          'original_exact_identity_benchmark':'preserved separately; not replaced',
          'no_screen_scores_read':True}
    (OUT/'analysis_plan.json').write_text(json.dumps(plan,indent=2));print(json.dumps(plan,indent=2),flush=True)
    sys.exit(0)
sys.path.insert(0, str(ROOT/'scripts'))
from product_execute_hiphop import load_task_relevance, rank_standardize_rows, compound_groups, select_best_conditions

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()

def log(*v):print(*v,flush=True)

protocol = {
 'question': 'Does target-specific strain ordering, beyond common assay coverage, support nominations and known-target retrieval?',
 'scope': 'Frozen 1175 targets, 3250 compounds, archived final task tables; no model or nomination changes.',
 'primary_randomization': 'A shared permutation of represented measured strain identities across all targets; preserves target-target profile relationships, fixed missing-strain tail and compound responses. Applied to final task profiles, not model retraining.',
 'n_mapping_permutations': 1000, 'seed': 20261011,
 'selection_replay': 'Reselect best condition and top 40 target-compound pairs in each direction; use archived stage-1 null_mean/null_sd as frozen score calibration for observed and randomized matrices.',
 'additional_null': '1000 shared permutations over all measured strains, matching archived exchangeability scope, for comparison.',
 'specificity': 'Decompose condition-level Spearman into across-target mean contribution plus target deviations. Descriptive comparisons of selected candidates are not independent validation.',
 'benchmark_identity': 'Exact full InChIKey overlap with archived Broad Repurposing Hub sample metadata; join drug annotations by pert_iname; all nonempty annotated targets intersected with the unchanged 1175-target universe.',
 'benchmark_primary': 'Direction-agnostic maximum of the two archived standardized directional scores; rank all 1175 targets per compound; macro-average known-target Recall@10, best-known-target reciprocal rank, and known-target rank percentile.',
 'benchmark_sensitivity': 'Report each direction separately and the uncorrected profile route, without changing primary endpoint after scores are read.',
 'benchmark_controls': 'Randomize target labels for annotated compounds; do not treat unannotated pairs as validated negatives. Existing four-pair panel retained separately.',
 'wet_experiments': False,
}
if represented_only:
    protocol.update(scope='Same frozen target and compound universes; score supplied task profiles only on represented measured strains.',
                    primary_randomization='Shared permutations of all represented measured strain identities; no unrepresented tail enters scoring.',
                    additional_null='Independent repeat of the same represented-strain permutation scope.',
                    selection_replay='Reselect best condition and top 40 pairs per direction using the corrected stage-1 null calibration.')
protocol['task_directory']=os.environ.get('YB_TASK_DIR',str(ROOT/'tasks/esm2_joint_tasks_dc_20260912'))
if os.environ.get('YB_STRAIN_STRATA'):
    protocol['primary_randomization']='Shared strain-label permutations within each published HIP/HOP assay stratum; preserves assay composition, target-target profile relationships and compound responses.'
    protocol['additional_null']='Independent repeat of the same within-assay permutation scope.'
    protocol['strain_strata_source']=os.environ['YB_STRAIN_STRATA']
    protocol['strain_strata_sha256']=sha(os.environ['YB_STRAIN_STRATA'])
external_protocol=os.environ.get('YB_EXTERNAL_PROTOCOL')
if external_protocol:
    protocol.update(json.loads(Path(external_protocol).read_text()))
    protocol['scope']='Frozen 1175 human target profiles; independent public chemical-genetic measurement universe'
    protocol['external_protocol_source']=external_protocol
    protocol['external_protocol_sha256']=sha(external_protocol)
(OUT/'analysis_plan.json').write_text(json.dumps(protocol,indent=2))
response_path=Path(os.environ.get('YB_RESPONSE_NPZ',ROOT/'raw/tier1_response/strain_response.npz'))
resp=np.load(response_path,allow_pickle=True)
orfs=resp['strain_orfs'].astype(str)
inchi=resp['compound_inchikeys'].astype(str)
keep=(~resp['is_vehicle'].astype(bool)) & (inchi!='')
conditions=np.flatnonzero(keep)
keys=inchi[keep]
groups=compound_groups(keys)
compounds=groups['ids']
pos={s:i for i,s in enumerate(orfs)}
task_dir=Path(os.environ.get('YB_TASK_DIR',ROOT/'tasks/esm2_joint_tasks_dc_20260912'))
files=sorted(task_dir.glob('yeast_task_*.tsv'))
targets=np.array([p.stem.replace('yeast_task_','') for p in files])
scope_path=os.environ.get('YB_TARGET_SCOPE')
if scope_path:
    target_scope=pd.read_csv(scope_path,sep='\t')
    assert target_scope.target_id.is_unique and set(target_scope.target_id)==set(targets)
    protocol.update(scope='Explicit target scope; all compounds retained',n_targets=len(targets),
                    target_scope_source=scope_path,target_scope_sha256=sha(scope_path),
                    benchmark_primary=f'Direction-agnostic maximum of both standardized directional scores; rank all {len(targets)} in-scope targets per compound.',
                    selection_replay='Full-screen maximum statistics across all pairs and both directions; nominee counts do not define the tested universe.')
    (OUT/'analysis_plan.json').write_text(json.dumps(protocol,indent=2))
task=[];supports=[];manifest=[]
for p in files:
    x,cov=load_task_relevance(p,pos,len(orfs));task.append(x)
    supports.append(np.ones(len(orfs),dtype=bool) if cov==len(orfs) else x>x.min())
    manifest.append({'path':str(p),'sha256':sha(p),'coverage':cov})
support=np.stack(supports)
assert np.all(support==support[0]), 'Coverage differs across targets; refine randomization.'
support=np.flatnonzero(support[0]);outside=np.setdiff1d(np.arange(len(orfs)),support)
task=np.stack(task);response=resp['z_score'][:,keep].T
if represented_only:
    task=task[:,support];response=response[:,support];orfs=orfs[support]
    support=np.arange(len(support));outside=np.array([],dtype=int)
task=rank_standardize_rows(task).astype(np.float32)
response=rank_standardize_rows(response).astype(np.float32)
log('Prepared ranks',task.shape,response.shape,'represented',len(support),'outside',len(outside))
np.savez_compressed(OUT/'rank_inputs.npz',task=task,response=response,targets=targets,compounds=compounds,orfs=orfs,conditions=conditions,condition_compounds=keys,support=support,group_order=groups['order'],group_starts=groups['starts'],group_counts=groups['counts'])

# The benchmark is assembled from external annotations before reading screen scores.
md=UP/'drug_ko_benchmark_v1/source_metadata'
sample=pd.read_csv(md/'repo-sample-annotation-20250818.txt',sep='\t',comment='!',dtype=str).fillna('')
drug=pd.read_csv(md/'repo-drug-annotation-20250818.txt',sep='\t',comment='!',dtype=str).fillna('')
matched=sample.loc[sample.InChIKey.isin(compounds),['InChIKey','pert_iname','pubchem_cid']].drop_duplicates().merge(drug,on='pert_iname',how='left').fillna('')
matched.to_csv(OUT/'benchmark_annotation_overlap.tsv',sep='\t',index=False)
rows=[]
for key,d in matched.groupby('InChIKey'):
    gold=sorted({g.strip() for v in d.target for g in re.split(r'[|;,]',v) if g.strip() in set(targets)})
    allgold=sorted({g.strip() for v in d.target for g in re.split(r'[|;,]',v) if g.strip()})
    for g in gold:
        rows.append({'inchikey':key,'target_id':g,'compound_names':'|'.join(sorted(set(d.pert_iname))), 'moa':'|'.join(sorted(set(d.moa))), 'all_annotated_targets':'|'.join(allgold), 'in_scope_gold_count':len(gold)})
gold=pd.DataFrame(rows)
if represented_only and not external_protocol:
    frozen=pd.read_csv(ROOT/'results_a6000/strengthening_20261011/benchmark_frozen_pairs.tsv',sep='\t',keep_default_na=False)
    if scope_path:
        frozen=frozen.loc[frozen.target_id.isin(targets)].reset_index(drop=True)
        frozen['in_scope_gold_count']=frozen.groupby('inchikey').target_id.transform('nunique')
    pd.testing.assert_frame_equal(gold,frozen)
if external_protocol and (OUT/'benchmark_frozen_pairs.tsv').exists():
    pd.testing.assert_frame_equal(gold,pd.read_csv(OUT/'benchmark_frozen_pairs.tsv',sep='\t',keep_default_na=False))
gold.to_csv(OUT/'benchmark_frozen_pairs.tsv',sep='\t',index=False)
log('Benchmark fixed',len(gold),'pairs',gold.inchikey.nunique(),'compounds',gold.target_id.nunique(),'targets')
if os.environ.get('YB_PREPARE_RANKS_ONLY')=='1':
    for p in [response_path,md/'repo-sample-annotation-20250818.txt',md/'repo-drug-annotation-20250818.txt',OUT/'benchmark_frozen_pairs.tsv']:
        manifest.append({'path':str(p),'sha256':sha(p)})
    (OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
    log('Ranks and benchmark frozen before screening evaluation');sys.exit(0)

# Exact saved screen matrices, retaining all compounds and hypotheses.
arrays={}
index=pd.MultiIndex.from_product([targets,compounds],names=['target_id','inchikey'])
for tag,folder in [('pos','results_esm2jointdc_screen_20260912'),('neg','results_esm2jointdc_neg_screen_20260912'),('raw','results_esm2joint_screen_20260912')]:
    if tag=='raw' and os.environ.get('YB_SKIP_RAW_SCREEN')=='1':continue
    p=UP/'product_pipeline/product/execute_hiphop'/folder/'exec_matrix.tsv'
    p=Path(os.environ.get('YB_'+tag.upper()+'_SCREEN',p))
    cols=['target_id','inchikey','spearman_rho','null_mean','null_sd','source_condition_index','emp_p']
    d=pd.read_csv(p,sep='\t',usecols=cols).set_index(['target_id','inchikey']).reindex(index)
    assert (scope_path or len(targets)==1175) and (external_protocol or len(compounds)==3250)
    assert len(d)==len(targets)*len(compounds) and not d.isna().any().any()
    for c in cols[2:]:arrays[tag+'_'+c]=d[c].to_numpy().reshape(len(targets),len(compounds)).astype(np.int32 if c=='source_condition_index' else np.float64)
    manifest.append({'path':str(p),'sha256':sha(p),'n_pairs':len(d)})
    log('Loaded',tag,len(d))
np.savez_compressed(OUT/'screen_arrays.npz',**arrays)
for p in [response_path, md/'repo-sample-annotation-20250818.txt',md/'repo-drug-annotation-20250818.txt',ROOT/'scripts/product_execute_hiphop.py']:
    manifest.append({'path':str(p),'sha256':sha(p)})
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
log('Preparation finished')
