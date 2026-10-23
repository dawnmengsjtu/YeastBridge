from __future__ import annotations
import json,sys,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path('/public/home/mengxl/dzy/YeastBridge_0912')
OUT=Path(os.environ.get('YB_ANALYSIS_DIR',ROOT/'results_a6000/strengthening_20261011'))
d=np.load(OUT/'rank_inputs.npz');s=dict(np.load(OUT/'screen_arrays.npz'))
targets=d['targets'];compounds=d['compounds'];T=len(targets)
tp={v:i for i,v in enumerate(targets)};cp={v:i for i,v in enumerate(compounds)}
gold=None if any(os.environ.get(k) for k in ['YB_CHEMBL_EVALUATE','YB_GPCR_PANEL_EVALUATE','YB_CANDIDATE_REPLICATION']) else pd.read_csv(OUT/'benchmark_frozen_pairs.tsv',sep='\t')
zp=(s['pos_spearman_rho']-s['pos_null_mean'])/s['pos_null_sd']
zn=(s['neg_spearman_rho']-s['neg_null_mean'])/s['neg_null_sd']
raw=(s['raw_spearman_rho']-s['raw_null_mean'])/s['raw_null_sd']
# Uncorrected archived data contain positive-direction maxima only: do not invent its negative-direction calibration.
methods={'final_bidirectional':np.maximum(zp,zn),'final_positive':zp,'final_negative':zn,'uncorrected_positive':raw}
if os.environ.get('YB_UNCORRECTED_NEGATIVE'):
    negative=np.load(os.environ['YB_UNCORRECTED_NEGATIVE'])
    assert np.array_equal(negative['targets'],targets) and np.array_equal(negative['compounds'],compounds)
    methods['uncorrected_negative']=(negative['spearman_rho']-negative['null_mean'])/negative['null_sd']
if os.environ.get('YB_CANDIDATE_REPLICATION'):
    from product_execute_hiphop import bh_fdr
    original=pd.read_csv(os.environ['YB_CANDIDATE_REPLICATION'],sep='\t')
    source_chem=pd.read_csv(os.environ['YB_CANDIDATE_CHEMISTRY'],sep='\t',keep_default_na=False).set_index('inchikey')
    external_chem=pd.read_csv(os.environ['YB_COMPOUND_TABLE'],sep='\t',keep_default_na=False)
    identity=os.environ.get('YB_REPLICATION_IDENTITY','exact');mapped=[]
    for row in original.itertuples():
        hits=[row.inchikey] if row.inchikey in cp else [];basis='exact'
        if not hits and identity=='parent':
            parent=source_chem.loc[row.inchikey,'parent_inchikey']
            if parent:hits=sorted(set(external_chem.loc[external_chem.parent_inchikey.eq(parent),'inchikey']) & set(cp))
            basis='parent'
        mapped.append(dict(direction=row.direction,target_id=row.target_id,original_inchikey=row.inchikey,external_inchikey=hits[0] if len(hits)==1 else '',identity_basis=basis if len(hits)==1 else 'ambiguous' if hits else 'unmapped',n_matches=len(hits)))
    frame=pd.DataFrame(mapped);frame.to_csv(OUT/'candidate_replication_mapping.tsv',sep='\t',index=False)
    available=frame[frame.external_inchikey.ne('')].drop_duplicates(['direction','target_id','external_inchikey']).copy();tests=[]
    if len(available):
        ti=np.array([tp[t] for t in available.target_id]);ci=np.array([cp[k] for k in available.external_inchikey]);C=len(compounds)
        # Biological sign mapping is fixed by the published readout definitions.
        percentiles={di:((T-rankdata(-score,axis=0))/(T-1),(C-rankdata(-score,axis=1))/(C-1)) for di,score in [('pos',zp),('neg',zn)]}
        by_compound=[np.flatnonzero(available.external_inchikey.to_numpy()==k) for k in sorted(available.external_inchikey.unique())]
        by_target=[np.flatnonzero(available.target_id.to_numpy()==k) for k in sorted(available.target_id.unique())]
        direction=available.direction.to_numpy()
        def values(ix,endpoint):
            v=np.empty(len(available))
            for di in ['pos','neg']:
                mask=direction==di;v[mask]=percentiles[di][endpoint][ix[mask],ci[mask]]
            return v
        def macro(v,groups):return float(np.mean([v[ix].mean() for ix in groups]))
        groups=[by_compound,by_target];names=['compound_macro_target_percentile','target_macro_compound_percentile']
        observed=[macro(values(ti,j),groups[j]) for j in range(2)]
        for j,name in enumerate(names):available[name]=values(ti,j)
        scope=pd.read_csv(ROOT.parent/'yeastbridge_re_mvp_stage_20260901/product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t').set_index('target_id').reindex(targets)
        blocks=[np.flatnonzero(scope.target_family.to_numpy()==g) for g in sorted(scope.target_family.unique())]
        rng=np.random.default_rng(20261022);B=10000
        for mode in (['shared_target_rows','shared_within_gpcr_channel_rows'] if len(blocks)>1 else ['shared_target_rows']):
            null=np.empty((B,2))
            for b in range(B):
                perm=rng.permutation(T) if mode=='shared_target_rows' else np.arange(T)
                if mode!='shared_target_rows':
                    for block in blocks:perm[block]=rng.permutation(block)
                for j in range(2):null[b,j]=macro(values(perm[ti],j),groups[j])
            for j,name in enumerate(names):tests.append(dict(identity=identity,randomization=mode,metric=name,n_pairs=len(available),n_compounds=len(by_compound),n_targets=len(by_target),observed=observed[j],null_mean=null[:,j].mean(),p_higher=(1+(null[:,j]>=observed[j]).sum())/(B+1)))
        result=pd.DataFrame(tests);result['q_bh_all_tests']=bh_fdr(result.p_higher.to_numpy());result.to_csv(OUT/'candidate_replication_tests.tsv',sep='\t',index=False)
    available.to_csv(OUT/'candidate_replication_pair_scores.tsv',sep='\t',index=False)
    report=dict(identity=identity,n_original_pairs=len(frame),n_mapped_original_pairs=int(frame.external_inchikey.ne('').sum()),n_available_pairs=len(available),n_available_compounds=available.external_inchikey.nunique(),n_ambiguous_pairs=int(frame.identity_basis.eq('ambiguous').sum()),interpretation='Exploratory cross-experiment profile replication; both sources use negative values for sensitivity/depletion, so algebraic directions are matched without a sign flip; not direct binding validation.')
    (OUT/'candidate_replication_summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);sys.exit(0)
if os.environ.get('YB_GPCR_PANEL_EVALUATE'):
    import subprocess
    from product_execute_hiphop import bh_fdr
    panel=pd.read_csv(os.environ['YB_GPCR_PANEL_EVALUATE'],sep='\t')
    panel=panel[panel.eligible].copy()
    chemistry=pd.read_csv(os.environ['YB_COMPOUND_TABLE'],sep='\t',keep_default_na=False)
    identity=os.environ.get('YB_GPCR_PANEL_IDENTITY','exact');links=[]
    for row in panel.drop_duplicates('catnr').itertuples():
        hits=[row.inchikey] if row.inchikey in cp else []
        basis='exact'
        if not hits and identity=='parent' and row.parent_inchikey:
            hits=sorted(set(chemistry.loc[chemistry.parent_inchikey.eq(row.parent_inchikey),'inchikey']) & set(cp));basis='parent'
        links.append(dict(catnr=row.catnr,source_inchikey=row.inchikey,screen_inchikey=hits[0] if len(hits)==1 else '',identity_basis=basis if len(hits)==1 else 'ambiguous' if hits else 'unmapped',n_screen_matches=len(hits)))
    links=pd.DataFrame(links);links.to_csv(OUT/'functional_identity_mapping.tsv',sep='\t',index=False)
    current=panel.merge(links[links.screen_inchikey.ne('')],on='catnr',validate='many_to_one')
    keys=sorted(current.catnr.unique());receptors=sorted(current.target_id.unique());n=len(keys)
    if n==0:raise ValueError('no uniquely mapped functional assay compounds')
    y=current.pivot(index='target_id',columns='catnr',values='activation_hit').reindex(index=receptors,columns=keys).to_numpy(dtype=bool)
    annotations=current.drop_duplicates('catnr').set_index('catnr').loc[keys]
    ci=[cp[k] for k in annotations.screen_inchikey];ti=[tp[t] for t in receptors]
    code='''import json,sys
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
mols=[Chem.MolFromSmiles(s) for s in json.load(sys.stdin)]
assert all(m is not None for m in mols)
print(json.dumps([MurckoScaffold.MurckoScaffoldSmiles(mol=m) or 'acyclic' for m in mols]))
'''
    scaffold=json.loads(subprocess.check_output([os.environ['YB_CHEMISTRY_PYTHON'],'-c',code],input=json.dumps(annotations.smiles.tolist()).encode()))
    blocks=[np.flatnonzero(np.asarray(scaffold)==v) for v in sorted(set(scaffold))]
    metrics=['AUROC','average_precision','recall_top10percent','enrichment_top10percent']
    valid=(y.sum(1)>0)&(y.sum(1)<n);rng=np.random.default_rng(20261021);B=10000
    perm=np.stack([rng.permutation(n) for _ in range(B)])
    boot=[np.concatenate([blocks[j] for j in rng.integers(len(blocks),size=len(blocks))]) for _ in range(B)]
    def evaluate(score,labels):
        labels=np.atleast_2d(labels).astype(float);nn=score.size;pos=labels.sum(1);den=pos*(nn-pos)
        ranks=rankdata(score,method='average');auc=np.divide(labels@ranks-pos*(pos+1)/2,den,out=np.full(len(labels),np.nan),where=den>0)
        order=np.argsort(-score,kind='stable');ends=np.r_[np.flatnonzero(score[order][1:]!=score[order][:-1]),nn-1]
        tpv=np.cumsum(labels[:,order],axis=1)[:,ends];gains=np.diff(np.c_[np.zeros(len(labels)),tpv],axis=1)
        ap=np.divide((gains*tpv/(ends+1)).sum(1),pos,out=np.full(len(labels),np.nan),where=pos>0)
        k=max(1,int(np.ceil(.1*nn)));cut=score[order[k-1]];above=score>cut;ties=score==cut
        weights=above.astype(float)+ties*(k-above.sum())/ties.sum();rec=np.divide(labels@weights,pos,out=np.full(len(labels),np.nan),where=pos>0)
        return np.column_stack([auc,ap,rec,rec/(k/nn)])
    tests=[];by_target=[];pair_rows=[];paired={}
    for method,fullscore in methods.items():
        score=fullscore[np.ix_(ti,ci)];point=np.array([evaluate(score[i],y[i])[0] for i in range(len(receptors))])
        for i,target in enumerate(receptors):
            by_target.append(dict(method=method,target_id=target,n_compounds=n,n_activation_hits=int(y[i].sum()),activation_prevalence=float(y[i].mean()),eligible_for_macro=bool(valid[i]),**dict(zip(metrics,point[i]))))
        for i,target in enumerate(receptors):
            for j,key in enumerate(keys):pair_rows.append(dict(method=method,target_id=target,catnr=key,inchikey=annotations.loc[key,'screen_inchikey'],score=float(score[i,j]),activation_hit=bool(y[i,j])))
        if not valid.any():continue
        null=np.mean([evaluate(score[i],y[i,perm]) for i in np.flatnonzero(valid)],axis=0)
        draws=[]
        for ix in boot:
            yy=y[valid][:,ix]
            if np.any((yy.sum(1)==0)|(yy.sum(1)==len(ix))):continue
            draws.append(np.mean([evaluate(score[i,ix],y[i,ix])[0] for i in np.flatnonzero(valid)],axis=0))
        distribution=np.asarray(draws);observed=point[valid].mean(0)
        paired[method]=(observed,distribution)
        for j,metric in enumerate(metrics):
            tests.append(dict(method=method,metric=metric,n_compounds=n,n_receptors=int(valid.sum()),n_scaffolds=len(blocks),observed=observed[j],bootstrap_low=np.quantile(distribution[:,j],.025) if len(draws)>1 else np.nan,bootstrap_high=np.quantile(distribution[:,j],.975) if len(draws)>1 else np.nan,n_valid_bootstraps=len(draws),null_mean=null[:,j].mean(),p_higher=(1+(null[:,j]>=observed[j]).sum())/(B+1)))
    frame=pd.DataFrame(tests)
    if len(frame):frame['q_bh_all_tests']=bh_fdr(frame.p_higher.to_numpy())
    frame.to_csv(OUT/'functional_panel_tests.tsv',sep='\t',index=False)
    pd.DataFrame(by_target).to_csv(OUT/'functional_panel_per_target.tsv',sep='\t',index=False)
    pd.DataFrame(pair_rows).to_csv(OUT/'functional_panel_pair_scores.tsv',sep='\t',index=False)
    if 'final_positive' in paired and 'uncorrected_positive' in paired:
        difference=paired['final_positive'][0]-paired['uncorrected_positive'][0]
        draws=paired['final_positive'][1]-paired['uncorrected_positive'][1]
        pd.DataFrame([dict(comparison='final_positive_minus_uncorrected_positive',metric=m,difference=difference[j],bootstrap_low=np.quantile(draws[:,j],.025),bootstrap_high=np.quantile(draws[:,j],.975),n_valid_bootstraps=len(draws)) for j,m in enumerate(metrics)]).to_csv(OUT/'functional_panel_correction_comparison.tsv',sep='\t',index=False)
    report=dict(identity=identity,n_compounds=n,n_receptors=len(receptors),n_eligible_receptors=int(valid.sum()),n_activation_hits=int(y.sum()),n_scaffolds=len(blocks),n_permutations=B,interpretation='Human GPCR activation in an engineered-yeast reporter; absent activation is not evidence of absent binding. Shared compound-label permutations and scaffold block resampling; all fixed directions retained.')
    report['bootstrap_scope']='Intervals condition on every eligible receptor retaining both classes in a scaffold draw; retained draw counts are reported. Primary inference uses all 10000 shared-label permutations.'
    report['correction_comparison_scope']='Exploratory matched positive-direction comparison added after initial panel outcomes; identical scaffold draws and all frozen models and identity rules retained.'
    (OUT/'functional_panel_summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);print(frame.to_string(index=False),flush=True)
    sys.exit(0)
if os.environ.get('YB_CHEMBL_TARGET_RANKING')=='1':
    from product_execute_hiphop import bh_fdr
    evidence=pd.read_csv(OUT/'experimental_pair_evidence.tsv',sep='\t')
    hub=set(pd.read_csv(OUT.parent/'parent_identity_benchmark/benchmark_frozen_pairs.tsv',sep='\t').inchikey)
    membership=pd.read_csv(os.environ['YB_PROFILE_FAMILIES'],sep='\t').set_index('target_id').profile_family
    plan={'question':'For a fixed human target, do yeast-based scores prioritize experimentally active compounds over measured weak compounds?',
          'scores':list(methods),'assays':['B','F'],'subsets':['all','outside_Hub_benchmark'],
          'unit':'target; each target requires both measured active and weak compounds',
          'interval':'10000 bootstrap draws of frozen profile-family groups; equal target averaging',
          'nulls':[f'shared permutation of target-score rows across all {T} in-scope targets','within-target labels randomized with shared compound random keys'],
          'multiplicity':'BH across all target-ranking tests and all previous 32 ChEMBL tests',
          'selection':'no scores, thresholds or test units selected from their outcomes'}
    (OUT/'target_ranking_plan.json').write_text(json.dumps(plan,indent=2))
    rng=np.random.default_rng(20261019);B=10000
    perm=np.stack([rng.permutation(T) for _ in range(B)])
    keys=rng.random((B,len(compounds)))
    tests=[];rows=[]
    for assay in ['B','F']:
        for subset in ['all','outside_Hub_benchmark']:
            cur=evidence[evidence.assay_type.eq(assay)]
            if subset!='all':cur=cur[~cur.inchikey.isin(hub)]
            units=[]
            for target,g in cur.groupby('target_id'):
                pos=[cp[k] for k in g[g.evidence_class.eq('active')].inchikey]
                neg=[cp[k] for k in g[g.evidence_class.eq('weak')].inchikey]
                if pos and neg:units.append((target,pos,neg))
            if not units:continue
            fam=np.array([membership.loc[t] for t,_,_ in units]);blocks=[np.flatnonzero(fam==f) for f in sorted(set(fam))]
            boot=[np.concatenate([blocks[j] for j in rng.integers(len(blocks),size=len(blocks))]) for _ in range(B)]
            for method,score in methods.items():
                values=[];row_null=np.zeros(B);label_null=np.zeros(B)
                for target,pos,neg in units:
                    ti=tp[target];a=score[ti,pos];b=score[ti,neg]
                    value=float(np.mean((a[:,None]>b[None,:])+.5*(a[:,None]==b[None,:])));values.append(value)
                    # Evaluate every shuffled target against the unchanged measured compound panel.
                    all_a=score[:,pos,None];all_b=score[:,None,neg]
                    all_auc=((all_a>all_b)+.5*(all_a==all_b)).mean((1,2))
                    row_null+=all_auc[perm[:,ti]]
                    rank=rankdata(score[ti,pos+neg],method='average');ix=np.argsort(keys[:,pos+neg],axis=1)[:,:len(pos)]
                    label_null+=(rank[ix].sum(1)-len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg))
                    rows.append(dict(assay_type=assay,subset=subset,method=method,target_id=target,profile_family=str(membership.loc[target]),n_active=len(pos),n_weak=len(neg),AUROC=value))
                values=np.array(values);draw=np.array([values[ix].mean() for ix in boot]);obs=values.mean()
                for name,null in [('shared_target_rows',row_null),('within_target_compound_labels',label_null)]:
                    null/=len(units)
                    tests.append(dict(assay_type=assay,subset=subset,method=method,null=name,n_targets=len(units),n_profile_families=len(blocks),n_pairs=sum(len(a)+len(b) for _,a,b in units),observed=obs,bootstrap_low=np.quantile(draw,.025) if len(blocks)>1 else np.nan,bootstrap_high=np.quantile(draw,.975) if len(blocks)>1 else np.nan,null_mean=null.mean(),p_higher=(1+(null>=obs).sum())/(B+1)))
    result=pd.DataFrame(tests);prior=pd.read_csv(OUT/'experimental_validation_tests.tsv',sep='\t')
    q=bh_fdr(np.r_[prior.p_higher,result.p_higher]);result['q_bh_combined_pharmacology']=q[len(prior):]
    prior['q_bh_combined_pharmacology']=q[:len(prior)]
    result.to_csv(OUT/'target_ranking_tests.tsv',sep='\t',index=False)
    prior.to_csv(OUT/'experimental_validation_combined_multiplicity.tsv',sep='\t',index=False)
    pd.DataFrame(rows).to_csv(OUT/'target_ranking_per_target.tsv',sep='\t',index=False)
    print(result.to_string(index=False),flush=True);sys.exit(0)
if os.environ.get('YB_CHEMBL_EVALUATE'):
    import subprocess
    from product_execute_hiphop import bh_fdr
    from scipy.stats import mannwhitneyu
    evidence=pd.read_csv(OUT/'experimental_pair_evidence.tsv',sep='\t')
    hub=set(pd.read_csv(OUT.parent/'parent_identity_benchmark/benchmark_frozen_pairs.tsv',sep='\t').inchikey)
    chemistry=pd.read_csv(OUT.parent/'parent_identity_benchmark/screen_compound_parent_identity.tsv',sep='\t').set_index('inchikey')
    keys=sorted(evidence.inchikey.unique())
    scaffold_code='''import json,sys
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
mols=[Chem.MolFromSmiles(s) for s in json.load(sys.stdin)]
assert all(m is not None for m in mols)
print(json.dumps([MurckoScaffold.MurckoScaffoldSmiles(mol=m) or 'acyclic' for m in mols]))
'''
    scaffolds=json.loads(subprocess.check_output([os.environ['YB_CHEMISTRY_PYTHON'],'-c',scaffold_code],input=json.dumps(chemistry.loc[keys,'smiles'].tolist()).encode()))
    scaffold_map=dict(zip(keys,scaffolds));rng=np.random.default_rng(20261017);B=10000
    rankmats={n:rankdata(-s,axis=0,method='average') for n,s in methods.items()}
    summaries=[];compound_rows=[];pair_rows=[]
    for assay_type in ['B','F']:
        for subset in ['all','outside_Hub_benchmark']:
            current=evidence[evidence.assay_type.eq(assay_type)].copy()
            if subset=='outside_Hub_benchmark':current=current[~current.inchikey.isin(hub)]
            active=current[current.evidence_class.eq('active')].sort_values(['inchikey','target_id'])
            if active.empty:continue
            gt=np.array([tp[t] for t in active.target_id]);gc=np.array([cp[c] for c in active.inchikey])
            group_keys=sorted(active.inchikey.unique());groups=[np.flatnonzero(active.inchikey.to_numpy()==k) for k in group_keys]
            starts=np.array([g[0] for g in groups]);counts=np.array([len(g) for g in groups])
            blocks=[np.flatnonzero(np.array([scaffold_map[k] for k in group_keys])==s) for s in sorted({scaffold_map[k] for k in group_keys})]
            bootstrap=[np.concatenate([blocks[j] for j in rng.integers(len(blocks),size=len(blocks))]) for _ in range(B)]
            nulls={n:np.empty(B) for n in methods}
            for b in range(B):
                perm=rng.permutation(T)
                for name,r in rankmats.items():nulls[name][b]=(np.add.reduceat((T-r[perm[gt],gc])/(T-1),starts)/counts).mean()
            for name,ranks in rankmats.items():
                values=np.add.reduceat((T-ranks[gt,gc])/(T-1),starts)/counts
                draws=np.array([values[ix].mean() for ix in bootstrap]);obs=values.mean()
                summaries.append({'assay_type':assay_type,'subset':subset,'method':name,'metric':'active_target_percentile',
                                  'n_compounds':len(groups),'n_pairs':len(active),'n_scaffolds':len(blocks),'observed':obs,
                                  'bootstrap_low':np.quantile(draws,.025),'bootstrap_high':np.quantile(draws,.975),
                                  'null_mean':nulls[name].mean(),'p_higher':(1+(nulls[name]>=obs).sum())/(B+1)})
                for key,v in zip(group_keys,values):compound_rows.append(dict(assay_type=assay_type,subset=subset,method=name,inchikey=key,metric='active_target_percentile',value=v,scaffold=scaffold_map[key]))
                if subset=='all':
                    for record,rk in zip(active.to_dict('records'),ranks[gt,gc]):pair_rows.append(dict(record,method=name,**{f'rank_among_{T}':rk}))
                auc=[];auc_keys=[];auc_null=np.zeros(B)
                for key,g in current.groupby('inchikey'):
                    pos=[tp[t] for t in g[g.evidence_class.eq('active')].target_id]
                    neg=[tp[t] for t in g[g.evidence_class.eq('weak')].target_id]
                    if not pos or not neg:continue
                    j=cp[key];v=mannwhitneyu(methods[name][pos,j],methods[name][neg,j]).statistic/(len(pos)*len(neg))
                    auc.append(v);auc_keys.append(key)
                    scores=methods[name][pos+neg,j];label_ranks=rankdata(scores,method='average');n_pos=len(pos)
                    label_rng=np.random.default_rng(20261018+cp[key])
                    for b in range(B):
                        ii=label_rng.choice(len(scores),n_pos,replace=False)
                        auc_null[b]+=(label_ranks[ii].sum()-n_pos*(n_pos+1)/2)/(n_pos*len(neg))
                    compound_rows.append(dict(assay_type=assay_type,subset=subset,method=name,inchikey=key,metric='measured_active_vs_weak_AUROC',value=v,scaffold=scaffold_map[key]))
                if auc:
                    values=np.array(auc);ablocks=[np.flatnonzero(np.array([scaffold_map[k] for k in auc_keys])==s) for s in sorted({scaffold_map[k] for k in auc_keys})]
                    draws=np.array([values[np.concatenate([ablocks[j] for j in rng.integers(len(ablocks),size=len(ablocks))])].mean() for _ in range(B)])
                    auc_null/=len(auc)
                    summaries.append(dict(assay_type=assay_type,subset=subset,method=name,metric='measured_active_vs_weak_AUROC',n_compounds=len(auc),n_scaffolds=len(ablocks),observed=values.mean(),bootstrap_low=np.quantile(draws,.025) if len(ablocks)>1 else np.nan,bootstrap_high=np.quantile(draws,.975) if len(ablocks)>1 else np.nan,null_mean=auc_null.mean(),p_higher=(1+(auc_null>=values.mean()).sum())/(B+1),leave_one_compound_out_min=(values.sum()-values.max())/(len(values)-1) if len(values)>1 else np.nan,leave_one_compound_out_max=(values.sum()-values.min())/(len(values)-1) if len(values)>1 else np.nan))
    result=pd.DataFrame(summaries);mask=result.p_higher.notna();result.loc[mask,'q_bh_all_tests']=bh_fdr(result.loc[mask,'p_higher'].to_numpy())
    result.to_csv(OUT/'experimental_validation_tests.tsv',sep='\t',index=False)
    pd.DataFrame(compound_rows).to_csv(OUT/'experimental_validation_per_compound.tsv',sep='\t',index=False)
    pd.DataFrame(pair_rows).to_csv(OUT/'experimental_active_pair_ranks.tsv',sep='\t',index=False)
    print(result.to_string(index=False),flush=True);sys.exit(0)
if os.environ.get('YB_DIAGNOSTIC_SCORES'):
    extra=dict(np.load(os.environ['YB_DIAGNOSTIC_SCORES']))
    diagnostic_names=['rank_raw','top25','top50','top100','top200'] if os.environ.get('YB_SPARSE_DIAGNOSTICS')=='1' else ['rank_raw','rank_vehicle','value_raw','value_vehicle']
    if os.environ.get('YB_BACKGROUND_DIAGNOSTICS')=='1':diagnostic_names=['rank_raw','target_standardized','target_percentile','double_centered','compound_then_target']
    if os.environ.get('YB_DIAGNOSTIC_NAMES'):
        diagnostic_names=json.loads(os.environ['YB_DIAGNOSTIC_NAMES'])
        if not diagnostic_names or len(set(diagnostic_names))!=len(diagnostic_names):raise ValueError('empty or duplicate diagnostic methods')
    for name in diagnostic_names:
        methods[name+'_positive']=extra[name+'_pos']
        methods[name+'_negative']=extra[name+'_neg']
        methods[name+'_bidirectional']=np.maximum(extra[name+'_pos'],extra[name+'_neg'])
if os.environ.get('YB_PREVIOUS_SCREENS'):
    previous=dict(np.load(os.environ['YB_PREVIOUS_SCREENS']))
    methods['previous_bidirectional']=np.maximum(*[(previous[di+'_spearman_rho']-previous[di+'_null_mean'])/previous[di+'_null_sd'] for di in ['pos','neg']])
metric_names=['recall_at_10','recall_at_50','best_known_target_RR','mean_known_target_percentile']
if os.environ.get('YB_WITHIN_FAMILY'):
    sys.path.insert(0,str(ROOT/'scripts'))
    from product_execute_hiphop import bh_fdr
    membership=pd.read_csv(os.environ['YB_WITHIN_FAMILY'],sep='\t').set_index('target_id').reindex(targets)
    blocks=[np.flatnonzero(membership.profile_family.to_numpy()==f) for f in sorted(membership.profile_family.unique())]
    units=[]
    for key,g in gold.groupby('inchikey',sort=True):
        known={tp[t] for t in g.target_id}
        for ix in blocks:
            hit=sorted(known & set(ix))
            if hit and 0<len(hit)<len(ix):
                for t in hit:units.append({'inchikey':key,'target_id':targets[t],'ti':t,'ci':cp[key],'family':membership.profile_family.iloc[t],'family_size':len(ix)})
    eligible=pd.DataFrame(units).sort_values(['inchikey','family','target_id']).reset_index(drop=True)
    if eligible.empty:raise ValueError('no nontrivial within-family known-target comparisons')
    eligible.to_csv(OUT/'within_family_eligible_pairs.tsv',sep='\t',index=False)
    gt=eligible.ti.to_numpy();gc=eligible.ci.to_numpy();sizes=eligible.family_size.to_numpy()
    groups=[np.flatnonzero(eligible.inchikey.to_numpy()==c) for c in sorted(eligible.inchikey.unique())]
    starts=np.array([g[0] for g in groups]);counts=np.array([len(g) for g in groups])
    def metric(rr):
        percentile=(sizes-rr)/(sizes-1)
        return np.array([(np.add.reduceat(percentile,starts)/counts).mean(),(1/np.minimum.reduceat(rr,starts)).mean()])
    family_ranks={}
    for name,score in methods.items():
        r=np.empty_like(score)
        for ix in blocks:r[ix]=rankdata(-score[ix],axis=0,method='average')
        family_ranks[name]=r
    rng=np.random.default_rng(20261016);B=10000;nulls={n:np.empty((B,2)) for n in family_ranks}
    for b in range(B):
        perm=np.arange(T)
        for ix in blocks:
            if len(ix)>1:perm[ix]=rng.permutation(ix)
        for name,r in family_ranks.items():nulls[name][b]=metric(r[perm[gt],gc])
    results=[];pair_results=[]
    for name,r in family_ranks.items():
        rr=r[gt,gc];obs=metric(rr);null=nulls[name]
        for k,m in enumerate(['mean_within_family_percentile','best_within_family_known_target_RR']):
            results.append({'method':name,'metric':m,'n_compounds':len(groups),'n_pairs':len(eligible),
                            'n_families':eligible.family.nunique(),'observed':obs[k],'null_mean':null[:,k].mean(),
                            'p_higher':float((1+(null[:,k]>=obs[k]).sum())/(B+1))})
        for row,v in zip(eligible.to_dict('records'),rr):pair_results.append(dict(row,method=name,within_family_rank=v))
    result=pd.DataFrame(results);result['q_bh_all']=bh_fdr(result.p_higher.to_numpy())
    result.to_csv(OUT/'within_family_tests.tsv',sep='\t',index=False)
    pd.DataFrame(pair_results).to_csv(OUT/'within_family_pair_ranks.tsv',sep='\t',index=False)
    print(result.to_string(index=False),flush=True);sys.exit(0)
rows=[];pairrows=[];per={};rankmats={}
for method,score in methods.items():
    method_pos,method_neg=zp,zn
    if method=='previous_bidirectional':
        method_pos,method_neg=[(previous[di+'_spearman_rho']-previous[di+'_null_mean'])/previous[di+'_null_sd'] for di in ['pos','neg']]
    if os.environ.get('YB_DIAGNOSTIC_SCORES'):
        prefix=method.rsplit('_',1)[0]
        if prefix+'_pos' in extra:method_pos,method_neg=extra[prefix+'_pos'],extra[prefix+'_neg']
    ranks=rankdata(-score,axis=0,method='average');rankmats[method]=ranks
    for key,g in gold.groupby('inchikey',sort=True):
        j=cp[key];ix=np.array([tp[t] for t in g.target_id]);rr=ranks[ix,j]
        r={'method':method,'inchikey':key,'compound_names':g.compound_names.iloc[0],'n_known_targets':len(ix),'best_known_target_rank':float(rr.min()),'recall_at_10':float((rr<=10).mean()),'recall_at_50':float((rr<=50).mean()),'best_known_target_RR':float(1/rr.min()),'mean_known_target_percentile':float(((T-rr)/(T-1)).mean()),'at_least_one_known_target_in_top10':bool(rr.min()<=10)}
        rows.append(r)
        for gi,rk in zip(ix,rr):
            direction='positive' if method.endswith('_positive') or (not method.endswith('_negative') and method_pos[gi,j]>=method_neg[gi,j]) else 'negative'
            pairrows.append({'method':method,'inchikey':key,'compound_names':g.compound_names.iloc[0],'target_id':targets[gi],f'rank_among_{T}':float(rk),'score':float(score[gi,j]),'preferred_direction':direction})
    per[method]=pd.DataFrame([r for r in rows if r['method']==method]).sort_values('inchikey')
pd.DataFrame(rows).to_csv(OUT/'known_target_per_compound.tsv',sep='\t',index=False)
pd.DataFrame(pairrows).to_csv(OUT/'known_target_pair_ranks.tsv',sep='\t',index=False)

rng=np.random.default_rng(20261012);B=10000;boot=rng.integers(0,len(per['final_bidirectional']),size=(B,len(per['final_bidirectional'])))
summary=[]
for method,frame in per.items():
    for metric in metric_names:
        v=frame[metric].to_numpy();draw=v[boot].mean(1)
        summary.append({'method':method,'metric':metric,'macro_mean':float(v.mean()),'bootstrap_95_low':float(np.quantile(draw,.025)),'bootstrap_95_high':float(np.quantile(draw,.975))})
pd.DataFrame(summary).to_csv(OUT/'known_target_macro_metrics.tsv',sep='\t',index=False)
if os.environ.get('YB_EXTERNAL_PROTOCOL') or os.environ.get('YB_TARGET_SCOPE'):
    import subprocess
    chemistry=pd.read_csv(os.environ['YB_COMPOUND_TABLE'],sep='\t').set_index('inchikey')
    keys=per['final_bidirectional'].inchikey.tolist()
    code='''import json,sys
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
mols=[Chem.MolFromSmiles(s) for s in json.load(sys.stdin)]
assert all(m is not None for m in mols)
print(json.dumps([MurckoScaffold.MurckoScaffoldSmiles(mol=m) or 'acyclic' for m in mols]))
'''
    scaffolds=json.loads(subprocess.check_output([os.environ['YB_CHEMISTRY_PYTHON'],'-c',code],input=json.dumps(chemistry.loc[keys,'smiles'].tolist()).encode()))
    blocks=[np.flatnonzero(np.asarray(scaffolds)==sc) for sc in sorted(set(scaffolds))]
    srng=np.random.default_rng(20261020)
    draws=[np.concatenate([blocks[j] for j in srng.integers(len(blocks),size=len(blocks))]) for _ in range(B)]
    intervals=[]
    values={name:frame[metric_names].to_numpy() for name,frame in per.items()}
    values['final_positive_minus_uncorrected_positive']=values['final_positive']-values['uncorrected_positive']
    if 'uncorrected_negative' in values:
        values['final_negative_minus_uncorrected_negative']=values['final_negative']-values['uncorrected_negative']
    for name,v in values.items():
        distribution=np.array([v[ix].mean(0) for ix in draws])
        for i,metric in enumerate(metric_names):intervals.append(dict(method=name,metric=metric,observed=float(v[:,i].mean()),n_compounds=len(keys),n_scaffolds=len(blocks),bootstrap_low=np.quantile(distribution[:,i],.025),bootstrap_high=np.quantile(distribution[:,i],.975)))
    pd.DataFrame(intervals).to_csv(OUT/'known_target_scaffold_intervals.tsv',sep='\t',index=False)
    pd.DataFrame({'inchikey':keys,'scaffold':scaffolds}).to_csv(OUT/'known_target_scaffolds.tsv',sep='\t',index=False)

# Paired positive-direction correction comparison; use the identical eligible compound set.
comparisons=[]
for metric in metric_names:
    diff=per['final_positive'][metric].to_numpy()-per['uncorrected_positive'][metric].to_numpy()
    draw=diff[boot].mean(1)
    comparisons.append({'comparison':'final_positive_minus_uncorrected_positive','metric':metric,'mean_difference':float(diff.mean()),'bootstrap_95_low':float(np.quantile(draw,.025)),'bootstrap_95_high':float(np.quantile(draw,.975))})
pd.DataFrame(comparisons).to_csv(OUT/'known_target_correction_comparison.tsv',sep='\t',index=False)
if 'previous_bidirectional' in per:
    changes=[]
    for metric in metric_names:
        diff=per['final_bidirectional'][metric].to_numpy()-per['previous_bidirectional'][metric].to_numpy()
        draw=diff[boot].mean(1)
        changes.append({'metric':metric,'mean_difference':float(diff.mean()),'bootstrap_95_low':float(np.quantile(draw,.025)),'bootstrap_95_high':float(np.quantile(draw,.975))})
    pd.DataFrame(changes).to_csv(OUT/'known_target_previous_comparison.tsv',sep='\t',index=False)

# Randomized known-target labels preserve every compound's number of known targets.
null=np.zeros((B,len(metric_names)),dtype=float)
ranks=rankmats['final_bidirectional']
for key,g in gold.groupby('inchikey',sort=True):
    j=cp[key];k=len(g)
    for b in range(B):
        rr=ranks[rng.choice(T,k,replace=False),j]
        null[b]+=[(rr<=10).mean(),(rr<=50).mean(),1/rr.min(),((T-rr)/(T-1)).mean()]
null/=gold.inchikey.nunique()
observed=per['final_bidirectional'][metric_names].mean().to_numpy()
nullout=[]
for k,metric in enumerate(metric_names):
    nullout.append({'metric':metric,'observed':float(observed[k]),'random_label_mean':float(null[:,k].mean()),'random_label_95_low':float(np.quantile(null[:,k],.025)),'random_label_95_high':float(np.quantile(null[:,k],.975)),'empirical_p_higher':float((1+(null[:,k]>=observed[k]).sum())/(B+1))})
pd.DataFrame(nullout).to_csv(OUT/'known_target_random_label_comparison.tsv',sep='\t',index=False)
if os.environ.get('YB_EXTERNAL_PROTOCOL') or os.environ.get('YB_TARGET_SCOPE'):
    # Same shared-label controls as strengthening_validate, without requiring discovery nominees.
    from product_execute_hiphop import bh_fdr
    ordered=gold.sort_values(['inchikey','target_id'])
    gt=np.array([tp[t] for t in ordered.target_id]);gc=np.array([cp[k] for k in ordered.inchikey])
    starts=np.r_[0,np.flatnonzero(ordered.inchikey.to_numpy()[1:]!=ordered.inchikey.to_numpy()[:-1])+1]
    counts=np.diff(np.r_[starts,len(ordered)])
    def calc(rr):
        return np.array([(np.add.reduceat(rr<=10,starts,dtype=float)/counts).mean(),(np.add.reduceat(rr<=50,starts,dtype=float)/counts).mean(),(1/np.minimum.reduceat(rr,starts)).mean(),(np.add.reduceat((T-rr)/(T-1),starts)/counts).mean()])
    scope=pd.read_csv(ROOT.parent/'yeastbridge_re_mvp_stage_20260901/product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv',sep='\t').set_index('target_id').reindex(targets)
    strata=[np.flatnonzero(scope.target_family.to_numpy()==g) for g in sorted(scope.target_family.unique())]
    drng=np.random.default_rng(20261013);tests=[]
    for mode in (['shared_target_permutation','shared_within_gpcr_channel_permutation'] if len(strata)>1 else ['shared_target_permutation']):
        nulls={name:np.empty((B,4)) for name in rankmats}
        for b in range(B):
            if mode=='shared_target_permutation':perm=drng.permutation(T)
            else:
                perm=np.arange(T)
                for ix in strata:perm[ix]=drng.permutation(ix)
            for name,rmat in rankmats.items():nulls[name][b]=calc(rmat[perm[gt],gc])
        for name,rmat in rankmats.items():
            obs=calc(rmat[gt,gc]);null=nulls[name];pv=(1+(null>=obs).sum(0))/(B+1)
            for i,metric in enumerate(metric_names):tests.append(dict(randomization=mode,method=name,metric=metric,observed=obs[i],null_mean=null[:,i].mean(),p_higher=pv[i]))
    frame=pd.DataFrame(tests);frame['q_bh_all_32' if len(tests)==32 else 'q_bh_all_tests']=bh_fdr(frame.p_higher.to_numpy())
    frame.to_csv(OUT/'known_target_direction_sensitivity.tsv',sep='\t',index=False)

# Original four-pair panel, evaluated with the same current target-ranking endpoint.
pairs=pd.read_csv(ROOT/'panels/positive_control_panel_20260911/pair_allowlist.tsv',sep='\t')
if os.environ.get('YB_TARGET_SCOPE'):pairs=pairs.loc[pairs.target_id.isin(targets)]
old=[]
for r in pairs.itertuples():
    for method in methods:
        available=r.target_id in tp and r.inchikey in cp
        old.append({'target_id':r.target_id,'inchikey':r.inchikey,'method':method,'available_in_dataset':available,f'rank_among_{T}':float(rankmats[method][tp[r.target_id],cp[r.inchikey]]) if available else np.nan})
pd.DataFrame(old).to_csv(OUT/'original_positive_panel_ranks.tsv',sep='\t',index=False)
result={'n_compounds':int(gold.inchikey.nunique()),'n_known_pairs':len(gold),'n_known_targets':int(gold.target_id.nunique()),'ranking_universe':T,'primary_method':'final_bidirectional','n_randomizations':B,'n_bootstraps':B,'seed':20261012,'hit_compounds_at10':int(per['final_bidirectional'].at_least_one_known_target_in_top10.sum()),'random_label_comparison':nullout,'limitations':['Retrospective evaluation of frozen final task profiles; not newly held-out model development data.','Known target annotations are positive labels; unannotated pairs are not verified inactive.','The uncorrected archive has only a positive-direction screen; only matched positive-direction before/after comparison is reported.','No model fitting or parameter selection uses this new benchmark.']}
if os.environ.get('YB_DIAGNOSTIC_SCORES'):
    result['limitations'][-1] = 'Primary model weights and screening rule remain frozen. Diagnostic variants are selected using training-scaffold folds only, then evaluated on held-out scaffolds.'
if os.environ.get('YB_EXTERNAL_PROTOCOL'):
    result['external_measurement_protocol']=os.environ['YB_EXTERNAL_PROTOCOL']
    result['limitations'][0]='Independent public chemical-genetic measurements; model and scoring frozen before external evaluation. Previously evaluated compounds or target labels may recur and are identified separately.'
if os.environ.get('YB_UNCORRECTED_NEGATIVE'):
    result['limitations'][2]='Both uncorrected directions supplied in a disclosed post-hoc sign-completeness analysis; the prespecified primary benchmark remains separate.'
(OUT/'known_target_benchmark_summary.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2),flush=True)

if os.environ.get('YB_DIAGNOSTIC_SCORES'):
    from sklearn.model_selection import GroupKFold
    candidates=[name+'_bidirectional' for name in diagnostic_names]
    compounds_order=per[candidates[0]].inchikey.to_numpy()
    chemistry=pd.read_csv('/public/home/mengxl/dzy/yeastbridge_vs/data/external/paired_compounds/hiphop_lee2014_e_mtab_2391/processed/compounds.tsv.gz',sep='\t')
    chemistry=chemistry[chemistry.inchikey.isin(compounds_order)].drop_duplicates(['inchikey','smiles'])
    if chemistry.inchikey.duplicated().any():raise ValueError('ambiguous benchmark structure')
    smiles=chemistry.set_index('inchikey').smiles.reindex(compounds_order)
    import subprocess
    chemistry_python=os.environ.get('YB_CHEMISTRY_PYTHON',sys.executable)
    scaffold_code='''import json,sys
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
mols=[Chem.MolFromSmiles(s) for s in json.load(sys.stdin)]
assert all(m is not None for m in mols), 'invalid benchmark SMILES'
print(json.dumps([MurckoScaffold.MurckoScaffoldSmiles(mol=m) or 'acyclic' for m in mols]))
'''
    scaffolds=json.loads(subprocess.check_output([chemistry_python,'-c',scaffold_code],input=json.dumps(smiles.tolist()).encode()))
    tensor=np.stack([per[name][metric_names].to_numpy() for name in candidates])
    held=np.empty_like(tensor[0]);fold_of=np.empty(len(compounds_order),dtype=int);choices=[]
    for fold,(train,test) in enumerate(GroupKFold(n_splits=5).split(compounds_order,groups=scaffolds),1):
        training_rr=tensor[:,train,2].mean(1);chosen=int(np.argmax(training_rr))
        held[test]=tensor[chosen,test];fold_of[test]=fold
        choices.append({'fold':fold,'n_train':len(train),'n_test':len(test),'selected_method':candidates[chosen],
                        'training_best_known_target_RR':float(training_rr[chosen])})
    pd.DataFrame(choices).to_csv(OUT/'crossfit_method_choices.tsv',sep='\t',index=False)
    held_frame=pd.DataFrame(held,columns=metric_names)
    held_frame.insert(0,'inchikey',compounds_order);held_frame['scaffold']=scaffolds;held_frame['fold']=fold_of
    held_frame.to_csv(OUT/'crossfit_per_compound.tsv',sep='\t',index=False)
    # Resample chemical scaffold groups, keeping related compounds together.
    group_indices=[np.flatnonzero(np.asarray(scaffolds)==g) for g in sorted(set(scaffolds))]
    rng=np.random.default_rng(20261014);draws=np.empty((10000,len(metric_names)))
    diff=held-tensor[0]
    for b in range(len(draws)):
        ii=np.concatenate([group_indices[j] for j in rng.integers(len(group_indices),size=len(group_indices))])
        draws[b]=diff[ii].mean(0)
    comparison=[]
    for k,metric in enumerate(metric_names):
        comparison.append({'metric':metric,'baseline_mean':float(tensor[0,:,k].mean()),'crossfit_mean':float(held[:,k].mean()),
                           'mean_difference':float(diff[:,k].mean()),'scaffold_bootstrap_95_low':float(np.quantile(draws[:,k],.025)),
                           'scaffold_bootstrap_95_high':float(np.quantile(draws[:,k],.975))})
    pd.DataFrame(comparison).to_csv(OUT/'crossfit_comparison.tsv',sep='\t',index=False)
    print('Scaffold-held-out comparison:',json.dumps(comparison),flush=True)
