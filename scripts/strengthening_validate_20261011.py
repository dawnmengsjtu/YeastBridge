from __future__ import annotations
import json, sys, hashlib, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path('/public/home/mengxl/dzy/YeastBridge_0912')
UP=Path('/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901')
OUT=Path(os.environ.get('YB_ANALYSIS_DIR',ROOT/'results_a6000/strengthening_20261011'))
TASK_DIR=Path(os.environ.get('YB_TASK_DIR',ROOT/'tasks/esm2_joint_tasks_dc_20260912'))
sys.path.insert(0,str(ROOT/'scripts'))
from product_execute_hiphop import bh_fdr, load_task_relevance
RESPONSE_PATH=Path(os.environ.get('YB_RESPONSE_NPZ',ROOT/'raw/tier1_response/strain_response.npz'))
if os.environ.get('YB_EVIDENCE_OVERVIEW')=='1':
    external=ROOT/'results_a6000/piotrowski2017_external';functional=ROOT/'results_a6000/published_gpcr_functional_2025'
    OUT.mkdir(parents=True,exist_ok=True);families={};sources=[]
    def combine(name,paths,expected_rows):
        pieces=[]
        for p in paths:
            f=pd.read_csv(p,sep='\t');f.insert(0,'source',str(p.relative_to(ROOT)));pieces.append(f)
            sources.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(f),family=name))
        f=pd.concat(pieces,ignore_index=True);assert len(f)==expected_rows,(name,len(f),expected_rows)
        assert np.isfinite(f.p_higher).all() and f.p_higher.between(0,1).all()
        f['q_bh_expanded_family']=bh_fdr(f.p_higher.to_numpy());f['expanded_family_size']=len(f)
        f.to_csv(OUT/(name+'.tsv'),sep='\t',index=False);families[name]=dict(n_tests=len(f),minimum_p=float(f.p_higher.min()),minimum_expanded_q=float(f.q_bh_expanded_family.min()),n_q_le_005=int(f.q_bh_expanded_family.le(.05).sum()))
        return f
    known=[];pharmacology=[];signs=[]
    for variant in ['primary','full_measured','genetic_projection']:
        b=external/variant
        known.extend([b/'known_target_direction_sensitivity.tsv',b/'parent_identity_benchmark/known_target_direction_sensitivity.tsv'])
        known.extend(b/'subgroups'/identity/subset/'known_target_direction_sensitivity.tsv' for identity in ['exact','parent'] for subset in ['source_all_conditions_high_confidence','outside_original_library','inside_original_library'])
        pharmacology.extend([b/'chembl_validation/experimental_validation_tests.tsv',b/'chembl_validation/target_ranking_tests.tsv'])
        signs.extend(b/'sign_complete_sensitivity'/identity/'known_target_direction_sensitivity.tsv' for identity in ['exact','parent'])
    a=combine('external_known_target_tests',known,768);b=combine('external_pharmacology_tests',pharmacology,192)
    pooled=pd.concat([a,b],ignore_index=True);pooled['q_bh_all_external_960']=bh_fdr(pooled.p_higher.to_numpy());pooled.to_csv(OUT/'external_all_960_tests.tsv',sep='\t',index=False)
    combine('external_sign_completion_tests',signs,240)
    models=json.loads((functional/'evaluation_plan.json').read_text())['models']
    paths=[functional/model/identity/'functional_panel_tests.tsv' for model in models for identity in ['exact','parent']]
    combine('published_functional_tests',paths,224)
    replication=[]
    for model in ['formal1241','assay1241','assay5865']:
        for identity in ['exact','parent']:
            d=external/'candidate_replication'/model/identity;j=json.loads((d/'candidate_replication_summary.json').read_text())
            if j['n_available_pairs']:replication.append(d/'candidate_replication_tests.tsv')
    combine('candidate_replication_tests',replication,20)
    (OUT/'overview_manifest.json').write_text(json.dumps(dict(families=families,sources=sources,interpretation='Expanded multiplicity sensitivity includes every frozen model and identity/subset comparison. Preserve original primary endpoints and original within-analysis q values; no branch selected by significance.'),indent=2))
    print(json.dumps(families,indent=2),flush=True);sys.exit(0)
if os.environ.get('YB_TAG_AUDIT') == '1':
    from scipy.stats import spearmanr
    from product_execute_hiphop import rank_standardize_rows
    tag=dict(np.load(RESPONSE_PATH.parent/'tag_response.npz'))
    response=dict(np.load(RESPONSE_PATH,allow_pickle=True))
    assert np.array_equal(tag['array_names'],response['array_names'])
    data=dict(np.load(OUT/'rank_inputs.npz'))
    destination=OUT/'barcode_reliability';destination.mkdir(exist_ok=True)
    common=np.isfinite(tag['up']).all(1) & np.isfinite(tag['down']).all(1)
    scope=common & np.isin(tag['strain_orfs'],data['orfs'])
    condition=[]
    for name,mask in [('all_measured',common),('task_represented',scope)]:
        up=tag['up'][mask].T;down=tag['down'][mask].T
        valid=(np.ptp(up,axis=1)>0) & (np.ptp(down,axis=1)>0)
        cor=np.full(len(up),np.nan)
        cor[valid]=(rank_standardize_rows(up[valid])*rank_standardize_rows(down[valid])).mean(1)
        for j,v in enumerate(cor):condition.append({'scope':name,'array':tag['array_names'][j],
                'is_vehicle':bool(tag['is_vehicle'][j]),'n_strains':int(mask.sum()),'uptag_downtag_spearman':float(v)})
    pd.DataFrame(condition).to_csv(destination/'all_array_concordance.tsv',sep='\t',index=False)
    strain=[]
    for i in np.flatnonzero(common):
        for group,cols in [('vehicle',tag['is_vehicle']),('compound',~tag['is_vehicle'])]:
            strain.append({'orf':tag['strain_orfs'][i],'group':group,'n_arrays':int(cols.sum()),
                           'up_vehicle_intensity':float(tag['up_vehicle_intensity'][i]),
                           'down_vehicle_intensity':float(tag['down_vehicle_intensity'][i]),
                           'uptag_downtag_spearman':float(spearmanr(tag['up'][i,cols],tag['down'][i,cols]).statistic)})
    pd.DataFrame(strain).to_csv(destination/'all_strain_concordance.tsv',sep='\t',index=False)
    candidates=pd.read_csv(OUT/'candidate_score_decomposition.tsv',sep='\t')
    screen=dict(np.load(OUT/'screen_arrays.npz'));tp={t:i for i,t in enumerate(data['targets'])};cp={c:i for i,c in enumerate(data['compounds'])}
    gene=tag['strain_orfs'][scope];records=[]
    for row in candidates.itertuples():
        frame=pd.read_csv(TASK_DIR/('yeast_task_'+row.target_id+'.tsv'),sep='\t').set_index('yeast_gene')
        x=frame.ensemble_score.reindex(gene).to_numpy();j=int(screen[row.direction+'_source_condition_index'][tp[row.target_id],cp[row.inchikey]])
        sign=1 if row.direction=='pos' else -1
        rec={'target_id':row.target_id,'inchikey':row.inchikey,'direction':row.direction,'array':tag['array_names'][j],
             'n_strains_with_both_tags':len(gene),'selection':'current combined-tag nominees and source condition, unchanged'}
        for kind in ['up','down']:rec[kind+'_rho']=float(spearmanr(x,sign*tag[kind][scope,j]).statistic)
        records.append(rec)
    pd.DataFrame(records).to_csv(destination/'candidate_barcode_concordance.tsv',sep='\t',index=False)
    report={'n_strains_both_tags':int(common.sum()),'n_represented_both_tags':int(scope.sum()),
            'n_candidates':len(records),'both_directional_rho_positive':sum(r['up_rho']>0 and r['down_rho']>0 for r in records),
            'same_condition_both_tags':True,'interpretation':'Technical barcode concordance; nominees were selected using their mean, so this is not independent biological or target validation.'}
    (destination/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True);sys.exit(0)
if os.environ.get('YB_RESPONSE_AUDIT') == '1':
    from scipy.stats import spearmanr
    OUT.mkdir(parents=True,exist_ok=True)
    current=dict(np.load(RESPONSE_PATH,allow_pickle=True))
    previous=dict(np.load(ROOT/'raw/tier1_response/strain_response.npz',allow_pickle=True))
    base=Path('/public/home/mengxl/dzy/yeastbridge_vs/data/external/paired_compounds/hiphop_lee2014_e_mtab_2391')
    metadata=pd.read_csv(base/'raw/E-MTAB-2391.sdrf.txt',sep='\t',dtype=str).fillna('').apply(lambda c:c.str.strip())
    old_columns={name:i for i,name in enumerate(previous['array_names'])}
    columns={name:i for i,name in enumerate(current['array_names'])}
    screens=pd.read_csv(UP/'yeastbridge_re_mvp/assets/drug_ko_benchmark_v1/eligible_screen_manifest.tsv',sep='\t').set_index('SGTCid')
    archive=UP/'drug_ko_benchmark_v1/official_fd_recovery_20260902/wayback_files'
    records=[];pairs=[];unmatched=[];sources=[]
    for path in sorted(archive.glob('SGTC_*_*.txt')):
        sid,mode=path.stem.rsplit('_',1)
        if sid not in screens.index:
            unmatched.append({'path':str(path),'reason':'no frozen screen metadata'});continue
        item=screens.loc[sid]
        concentration=float(item.conc)/(1000 if item.unit=='nM' else 1)
        dose=pd.to_numeric(metadata['Factor Value[dose]'],errors='coerce')
        dose=dose/np.where(metadata['Unit[concentration unit]'].eq('nanomolar'),1000,1)
        matching=metadata[metadata['Factor Value[compound]'].str.casefold().eq(str(item.condition).strip().casefold()) & np.isclose(dose,concentration,rtol=.02,atol=.001)]
        if len(matching)!=1:
            unmatched.append({'path':str(path),'reason':f'expected one exact-name/dose CEL, found {len(matching)}'});continue
        array=matching.iloc[0]['Array Data File']
        author=pd.read_csv(path,sep='\t');score=[c for c in author if c.endswith('.score')]
        if len(score)!=1 or author.orf.duplicated().any():
            raise ValueError(f'ambiguous official profile: {path}')
        author=author.set_index('orf');score=score[0]
        x=pd.DataFrame({'official_FD':author[score]})
        for name,data,cols in [('original_ORF_probes',previous,old_columns),('corrected_barcodes',current,columns)]:
            x[name]=pd.Series(-data['z_score'][:,cols[array]],index=data['strain_orfs']).reindex(x.index)
        x=x.dropna()
        for name in ['original_ORF_probes','corrected_barcodes']:
            high=min(50,len(x));a=set(x.nlargest(high,'official_FD').index);b=set(x.nlargest(high,name).index)
            records.append({'screen':sid,'condition':item.condition,'mode':mode,'array':array,'method':name,
                            'n_common_strains':len(x),'spearman_vs_official':float(spearmanr(x.official_FD,x[name]).statistic),
                            'top50_overlap':len(a & b)})
        ranked=x.rank(ascending=False,method='average')
        for gene in x.nlargest(5,'official_FD').index:
            pairs.append({'screen':sid,'mode':mode,'orf':gene,'gene':author.loc[gene,'gene'],
                          'n_common_strains':len(x),'official_rank':ranked.loc[gene,'official_FD'],
                          'old_rank':ranked.loc[gene,'original_ORF_probes'],'corrected_rank':ranked.loc[gene,'corrected_barcodes']})
        sources.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    pd.DataFrame(records).to_csv(OUT/'official_response_concordance.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(OUT/'official_top5_recovery.tsv',sep='\t',index=False)
    # Independent classic yeast controls; all available concentrations, no best-condition selection.
    controls=[]
    for key,gene in [('RFHAOTPXVQNOHP-UHFFFAOYSA-N','YHR007C'),('BWMISRWJRUSYEX-SZKNIZGXSA-N','YGR175C')]:
        for j in np.flatnonzero(current['compound_inchikeys']==key):
            name=current['array_names'][j]
            for method,data,cols in [('original_ORF_probes',previous,old_columns),('corrected_barcodes',current,columns)]:
                i=list(data['strain_orfs']).index(gene);v=data['z_score'][:,cols[name]]
                controls.append({'inchikey':key,'orf':gene,'array':name,'dose':current['doses'][j],
                                 'method':method,'depletion_rank':float(rankdata(v)[i]),
                                 'n_strains':len(v),'barcode_contrast':float(v[i])})
    pd.DataFrame(controls).to_csv(OUT/'yeast_positive_control_recovery.tsv',sep='\t',index=False)
    report={'n_official_profiles_compared':len(records)//2,'n_screens':len({r['screen'] for r in records}),
            'unmatched_official_files':unmatched,'official_sources':sources,
            'response_path':str(RESPONSE_PATH),'response_sha256':hashlib.sha256(RESPONSE_PATH.read_bytes()).hexdigest(),
            'interpretation':'Assay reconstruction checks against published profiles, not new human-target validation.'}
    (OUT/'response_audit_summary.json').write_text(json.dumps(report,indent=2))
    print(pd.DataFrame(records).groupby('method')[['spearman_vs_official','top50_overlap']].median().to_string(),flush=True)
    print(json.dumps(report,indent=2),flush=True)
    sys.exit(0)
d=dict(np.load(OUT/'rank_inputs.npz'));s=dict(np.load(OUT/'screen_arrays.npz'))
targets=d['targets'];compounds=d['compounds'];T=len(targets)
tp={v:i for i,v in enumerate(targets)};cp={v:i for i,v in enumerate(compounds)}
gold=pd.read_csv(OUT/'benchmark_frozen_pairs.tsv',sep='\t').sort_values(['inchikey','target_id'])
zp=(s['pos_spearman_rho']-s['pos_null_mean'])/s['pos_null_sd'];zn=(s['neg_spearman_rho']-s['neg_null_mean'])/s['neg_null_sd']
ranks=rankdata(-np.maximum(zp,zn),axis=0)
gt=np.array([tp[x] for x in gold.target_id]);gc=np.array([cp[x] for x in gold.inchikey])
starts=np.r_[0,np.flatnonzero(gold.inchikey.to_numpy()[1:]!=gold.inchikey.to_numpy()[:-1])+1]
counts=np.diff(np.r_[starts,len(gold)])
metrics=['recall_at_10','recall_at_50','best_known_target_RR','mean_known_target_percentile']
def calc(rr):
    return np.array([(np.add.reduceat(rr<=10,starts,dtype=float)/counts).mean(),(np.add.reduceat(rr<=50,starts,dtype=float)/counts).mean(),(1/np.minimum.reduceat(rr,starts)).mean(),(np.add.reduceat((T-rr)/(T-1),starts)/counts).mean()])
observed=calc(ranks[gt,gc]);rng=np.random.default_rng(20261013);B=10000
scope=pd.read_csv(UP/'product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t').set_index('target_id').reindex(targets)
strata=[np.flatnonzero(scope.target_family.to_numpy()==g) for g in sorted(scope.target_family.unique())]
rows=[]
for mode in ['shared_target_permutation','shared_within_gpcr_channel_permutation']:
    null=np.empty((B,4))
    for b in range(B):
        if mode=='shared_target_permutation':perm=rng.permutation(T)
        else:
            perm=np.arange(T)
            for ix in strata:perm[ix]=rng.permutation(ix)
        null[b]=calc(ranks[perm[gt],gc])
    pv=(1+(null>=observed).sum(0))/(B+1);qv=bh_fdr(pv)
    for i,m in enumerate(metrics):
        rows.append({'randomization':mode,'metric':m,'observed':observed[i],'null_mean':null[:,i].mean(),'null_95_low':np.quantile(null[:,i],.025),'null_95_high':np.quantile(null[:,i],.975),'p_higher':pv[i],'q_bh_four_metrics':qv[i]})
pd.DataFrame(rows).to_csv(OUT/'known_target_dependence_sensitivity.tsv',sep='\t',index=False)

# Reattach frozen task-family membership to all-target C1/C2 context.
famdir=Path(os.environ.get('YB_FAMILY_PANEL_DIR',ROOT/'panels/family_specificity_panel_20260913'))
membership=None
for path in famdir.iterdir():
    if path.suffix in ['.tsv','.csv'] and 'cluster' in path.name:
        a=pd.read_csv(path,sep='\t' if path.suffix=='.tsv' else ',')
        if 'target_id' in a.columns:
            membership=(path,a);break
if membership is None:
    # Reproduce the archived average-linkage family definition, without running its file-writing script.
    from scipy.cluster.hierarchy import fcluster,linkage
    from scipy.spatial.distance import squareform
    rec=json.loads((TASK_DIR/'export_record.json').read_text())
    tail=set(rec['tail_strains']);mat=[];orfs=None
    for t in targets:
        frame=pd.read_csv(TASK_DIR/('yeast_task_'+t+'.tsv'),sep='\t').set_index('yeast_gene')
        if orfs is None:orfs=[g for g in frame.index if g not in tail]
        mat.append(frame.ensemble_score.reindex(orfs).to_numpy())
    dist=1-np.corrcoef(np.asarray(mat));np.fill_diagonal(dist,0);dist=(dist+dist.T)/2
    lab=fcluster(linkage(squareform(dist,checks=False),method='average'),t=.1,criterion='distance')
    prev=json.loads((famdir/'build_record.json').read_text())
    assert len(set(lab))==prev['n_clusters_0.9']
    for fid,info in prev['families'].items():
        ids=[int(lab[tp[t]]) for t in info['nominated']]
        assert len(set(ids))==1 and int((lab==ids[0]).sum())==info['size']
    a=pd.DataFrame({'target_id':targets,'profile_family':['FAM%05d'%i for i in lab]})
    path=OUT/'reproduced_profile_family_membership.tsv';a.to_csv(path,sep='\t',index=False)
    membership=(path,a)
context=pd.read_csv(OUT/'C1_C2_all_target_scores.tsv',sep='\t')
context=context.drop(columns=['target_family','comparison_group'],errors='ignore')
context=context.merge(scope[['target_family']],left_on='target_id',right_index=True,how='left',validate='many_to_one')
if membership:
    path,a=membership
    cols=[c for c in a.columns if 'cluster' in c or 'family' in c]
    context=context.drop(columns=cols,errors='ignore')
    context=context.merge(a[['target_id']+cols],on='target_id',how='left',validate='many_to_one')
    if cols:
        c=cols[0]
        context['comparison_group']='other_targets'
        for label in ['C1','C2']:
            fams=set(context.loc[(context.compound==label)&context.nominated,c])
            context.loc[(context.compound==label)&context[c].isin(fams),'comparison_group']='same_profile_family_not_nominated'
        context.loc[context.nominated,'comparison_group']='nominated'
context.to_csv(OUT/'C1_C2_all_target_scores.tsv',sep='\t',index=False)
if 'comparison_group' in context:
    context.groupby(['compound','comparison_group']).agg(n=('target_id','size'),rho_min=('rho','min'),rho_median=('rho','median'),rho_max=('rho','max'),target_deviation_median=('target_deviation','median')).reset_index().to_csv(OUT/'C1_C2_context_summary.tsv',sep='\t',index=False)

# Summarize original selected pairs under both diagnostic nulls, retaining selection adjustment.
summary={}
for mode in ['represented_only','all_strains']:
    a=pd.read_csv(OUT/(mode+'_candidate_pvalues.tsv'),sep='\t')
    b=pd.read_csv(OUT/(mode+'_selection_replay.tsv'),sep='\t')
    summary[mode]={}
    for di in ['pos','neg']:
        x=a[a.direction==di];y=b[b.direction==di]
        summary[mode][di]={'fixed_candidate_pairs_nominal_p_le_005':int((x.p_fixed_pair<=.05).sum()),'fixed_candidate_p_min':float(x.p_fixed_pair.min()),'fixed_candidate_p_median':float(x.p_fixed_pair.median()),'screen_max_adjusted_p_min':float(x.p_screen_max_both_directions.min()),'screen_max_adjusted_pairs_p_le_005':int((x.p_screen_max_both_directions<=.05).sum()),'C1_mean_top40_count':float(y.C1_count.mean()),'C1_all40_frequency':float((y.C1_count==40).mean()),'C2_mean_top40_count':float(y.C2_count.mean()),'C2_at_least35_frequency':float((y.C2_count>=35).mean()),'largest_compound_mean_count':float(y.largest_compound_count.mean())}
(OUT/'mapping_control_summary.json').write_text(json.dumps(summary,indent=2))

# Independent numerical checks on individual pairs use SciPy Spearman rather than matrix multiplication.
from scipy.stats import spearmanr
resp=np.load(RESPONSE_PATH,allow_pickle=True)
pos={g:i for i,g in enumerate(resp['strain_orfs'].astype(str))}
checks=[]
for t,c in [('DRD5','SQFNCEKTEKNPCN-UHFFFAOYSA-N'),('ITPR3','CLQPEZWMMZFTBG-UHFFFAOYSA-N'),('DRD2','LNEPOXFFQSENCJ-UHFFFAOYSA-N')]:
    task,_=load_task_relevance(TASK_DIR/('yeast_task_'+t+'.tsv'),pos,len(pos))
    scoring_mask=np.isin(resp['strain_orfs'].astype(str),d['orfs'])
    ci=cp[c];ti=tp[t]
    for di,sign in [('pos',1),('neg',-1)]:
        col=int(s[di+'_source_condition_index'][ti,ci]);rho=float(spearmanr(task[scoring_mask],sign*resp['z_score'][scoring_mask,col]).statistic)
        archived=float(s[di+'_spearman_rho'][ti,ci]);checks.append({'target_id':t,'inchikey':c,'direction':di,'scipy_rho':rho,'archived_rho':archived,'absolute_error':abs(rho-archived)})
assert max(r['absolute_error'] for r in checks)<2e-7
(OUT/'independent_score_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps({'dependence_sensitivity':rows,'mapping_summary':summary,'family_membership_source':str(membership[0]) if membership else None,'independent_score_max_error':max(r['absolute_error'] for r in checks)},indent=2))

# Paired functional-query checks after an upstream input-contract correction.
if (OUT/'functional_v9/per_query_rr_basic_structure.tsv').exists():
    internal=pd.read_csv(OUT/'functional_v9/per_query_rr_basic_structure.tsv',sep='\t')
    external=pd.read_csv(OUT/'functional_ev1a/per_query_rr_ev1a.tsv',sep='\t')
    internal['unit']='Norman:'+internal.split_group.astype(str)
    external['unit']='EV1a:'+external.gene.astype(str)
    pooled=pd.concat([internal,external],ignore_index=True)
    frng=np.random.default_rng(20261012);permutations=100000;stats=[]
    for cohort,frame in [('Norman',internal),('EV1a',external),('pooled',pooled)]:
        for comparator in ['rr_human_only','rr_yeast_only','rr_comparator_best_single']:
            delta=frame.rr_joint_nested.to_numpy()-frame[comparator].to_numpy()
            for unit in ['query','biological_group']:
                values=delta if unit=='query' else pd.Series(delta).groupby(frame.unit.to_numpy()).sum().to_numpy()
                hits=0
                for start in range(0,permutations,1000):
                    signs=frng.choice([-1,1],size=(min(1000,permutations-start),len(values)))
                    hits+=int(((signs@values)>=delta.sum()).sum())
                stats.append({'cohort':cohort,'comparator':comparator,'randomization_unit':unit,
                              'n_queries':len(frame),'n_units':len(values),'delta_mrr':float(delta.mean()),
                              'n_permutations':permutations,'seed':20261012,
                              'one_sided_p_plus_one':float((hits+1)/(permutations+1))})
    pd.DataFrame(stats).to_csv(OUT/'functional_gain_sensitivity.tsv',sep='\t',index=False)

# Report every predefined direction, with joint multiplicity control, rather than selecting the best one.
raw_z=(s['raw_spearman_rho']-s['raw_null_mean'])/s['raw_null_sd']
rank_sets={'bidirectional':ranks,'positive':rankdata(-zp,axis=0),
           'negative':rankdata(-zn,axis=0),'uncorrected_positive':rankdata(-raw_z,axis=0)}
direction_rows=[];drng=np.random.default_rng(20261013)
for mode in ['shared_target_permutation','shared_within_gpcr_channel_permutation']:
    nulls={name:np.empty((B,4)) for name in rank_sets}
    for b in range(B):
        if mode=='shared_target_permutation':perm=drng.permutation(T)
        else:
            perm=np.arange(T)
            for ix in strata:perm[ix]=drng.permutation(ix)
        for name,rmat in rank_sets.items():nulls[name][b]=calc(rmat[perm[gt],gc])
    for name,rmat in rank_sets.items():
        obs=calc(rmat[gt,gc]);null=nulls[name];pv=(1+(null>=obs).sum(0))/(B+1)
        for i,metric in enumerate(metrics):
            direction_rows.append({'randomization':mode,'method':name,'metric':metric,
                                   'observed':obs[i],'null_mean':null[:,i].mean(),'p_higher':pv[i]})
direction_frame=pd.DataFrame(direction_rows)
direction_frame['q_bh_all_32']=bh_fdr(direction_frame.p_higher.to_numpy())
direction_frame.to_csv(OUT/'known_target_direction_sensitivity.tsv',sep='\t',index=False)
