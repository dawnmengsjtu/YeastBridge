from __future__ import annotations
import json, sys, time, os
from pathlib import Path
import numpy as np
import pandas as pd
import torch

ROOT=Path('/public/home/mengxl/dzy/YeastBridge_0912')
OUT=Path(os.environ.get('YB_ANALYSIS_DIR',ROOT/'results_a6000/strengthening_20261011'))
PANEL=Path(os.environ.get('YB_NOMINATION_DIR',ROOT/'panels/ej_dc_confirmation_panel_20260912'))
sys.path.insert(0,str(ROOT/'scripts'))
from product_execute_hiphop import compound_groups

torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32=False
torch.backends.cudnn.allow_tf32=False
d=dict(np.load(OUT/'rank_inputs.npz'));s=dict(np.load(OUT/'screen_arrays.npz'))
targets=d['targets'];compounds=d['compounds'];T=len(targets);C=len(compounds);S=d['task'].shape[1]
tp={v:i for i,v in enumerate(targets)};cp={v:i for i,v in enumerate(compounds)}
dev='cuda:0'
X=torch.tensor(d['task'],device=dev);Y=torch.tensor(d['response'].T,device=dev)
support=d['support'];outside=np.setdiff1d(np.arange(S),support)
strata_path=os.environ.get('YB_STRAIN_STRATA')
strata_labels=None
if strata_path:
    import hashlib
    frame=pd.read_csv(strata_path,sep='\t').set_index('orf')
    if not frame.index.is_unique:raise ValueError('duplicate strain annotation')
    labels=frame.assay.reindex(d['orfs'])
    if labels.isna().any():raise ValueError('missing assay stratum')
    strata_labels=labels.to_numpy()
    (OUT/'permutation_strata.json').write_text(json.dumps(dict(source=strata_path,sha256=hashlib.sha256(Path(strata_path).read_bytes()).hexdigest(),counts=labels.value_counts().to_dict(),rule='Shuffle only within published HIP/HOP assay; preserve assay membership and all original observed scores.'),indent=2))
def permute_strains(indices):
    if strata_labels is None:return rng.permutation(len(indices))
    labels=strata_labels[indices];perm=np.arange(len(indices))
    for label in sorted(set(labels)):
        take=np.flatnonzero(labels==label);perm[take]=rng.permutation(take)
    return perm
Xsub=X[:,support].contiguous();Ysub=Y[support].contiguous()
fixed=X[:,outside]@Y[outside]/S
groups=compound_groups(d['condition_compounds'])
order=torch.tensor(groups['order'],device=dev)
starts=torch.tensor(groups['starts'],device=dev)
multiple=[(i,int(a),int(n)) for i,(a,n) in enumerate(zip(groups['starts'],groups['counts'])) if n>1]

def reduce_scores(a):
    a=a[:,order]
    pos=a[:,starts].clone();neg=-a[:,starts].clone()
    for i,st,n in multiple:
        pos[:,i]=a[:,st:st+n].amax(1)
        neg[:,i]=(-a[:,st:st+n]).amax(1)
    return pos,neg

if os.environ.get('YB_BACKGROUND_DIAGNOSTICS') == '1':
    from scipy.stats import rankdata, norm, median_abs_deviation
    out=Path(os.environ['YB_DIAGNOSTIC_DIR']);out.mkdir(parents=True,exist_ok=True)
    plan={'question':'Do target background propensities and broad compound responses mask target ranking?',
          'variants':['rank_raw','target_standardized','target_percentile','double_centered','compound_then_target'],
          'normalization':'label-free, all frozen library scores; no bioactivity labels enter score transforms',
          'selection':'five-fold chemical-scaffold-held-out Hub evaluation; select using training best-known-target RR, retain baseline among choices',
          'test_metrics':['recall_at_10','recall_at_50','best_known_target_RR','mean_known_target_percentile'],
          'primary_screen':'unchanged until independently supported improvement',
          'motivation_reference':'10.1016/j.cell.2017.10.049; background-relative connectivity, adapted rather than identical to tau'}
    (out/'background_calibration_plan.json').write_text(json.dumps(plan,indent=2))
    def standardize(a,axis):
        center=np.median(a,axis=axis,keepdims=True)
        scale=median_abs_deviation(a,axis=axis,scale='normal',nan_policy='raise')
        if np.any(scale<=0) or np.any(~np.isfinite(scale)):raise ValueError('degenerate background scale')
        return (a-center)/np.expand_dims(scale,axis)
    result={}
    for di in ['pos','neg']:
        a=(s[di+'_spearman_rho']-s[di+'_null_mean'])/s[di+'_null_sd']
        variants={'rank_raw':a,'target_standardized':standardize(a,1),
                  'target_percentile':norm.ppf((rankdata(a,axis=1,method='average')-.5)/C),
                  'double_centered':a-a.mean(1,keepdims=True)-a.mean(0,keepdims=True)+a.mean(),
                  'compound_then_target':standardize(standardize(a,0),1)}
        for name,value in variants.items():
            assert np.isfinite(value).all();result[name+'_'+di]=value
    np.savez_compressed(out/'diagnostic_scores.npz',**result)
    print('Five background scoring variants saved; no screening nominees changed.',flush=True);sys.exit(0)

if os.environ.get('YB_SCORE_DIAGNOSTICS') == '1':
    from product_execute_hiphop import rank_standardize_rows
    from scipy.stats import median_abs_deviation
    out=Path(os.environ['YB_DIAGNOSTIC_DIR'])
    response=np.load(Path(os.environ.get('YB_RESPONSE_NPZ',ROOT/'raw/tier1_response/strain_response.npz')),allow_pickle=True)
    row_of={g:i for i,g in enumerate(response['strain_orfs'].astype(str))}
    ix=np.array([row_of[g] for g in d['orfs']]);z=response['z_score'][ix]
    vehicle=z[:,response['is_vehicle'].astype(bool)]
    center=np.median(vehicle,axis=1);scale=median_abs_deviation(vehicle,axis=1,scale='normal')
    if (np.any(~np.isfinite(scale)) or np.any(scale<=0)) and os.environ.get('YB_SPARSE_DIAGNOSTICS')!='1':
        raise ValueError('non-positive vehicle MAD; no automatic strain exclusion')
    pd.DataFrame({'yeast_gene':d['orfs'],'vehicle_median':center,'vehicle_mad':scale}).to_csv(out/'vehicle_calibration.tsv',sep='\t',index=False)
    if os.environ.get('YB_SPARSE_DIAGNOSTICS')=='1':
        Yvehicle=Y  # Sparse variants use the original responses; vehicle calibration is not an input.
    else:
        calibrated=(z[:,d['conditions']]-center[:,None])/scale[:,None]
        Yvehicle=torch.tensor(rank_standardize_rows(calibrated.T).T,dtype=torch.float32,device=dev)
    values=[]
    for target in targets:
        f=pd.read_csv(Path(os.environ.get('YB_TASK_DIR',ROOT/'tasks/esm2_joint_tasks_dc_20260912'))/('yeast_task_'+target+'.tsv'),sep='\t').set_index('yeast_gene')
        values.append(f.ensemble_score.reindex(d['orfs']).to_numpy())
    values=np.asarray(values);values-=values.mean(1,keepdims=True)
    if np.any(~np.isfinite(values)) or np.any(values.std(1)<=0):
        raise ValueError('non-finite or constant task values')
    values/=values.std(1,keepdims=True)
    Xvalue=torch.tensor(values,dtype=torch.float32,device=dev)
    seeds=np.random.default_rng(42).integers(0,np.iinfo(np.uint32).max,size=1000,dtype=np.uint32)
    diagnostic={};summary={}
    variants=[('rank_raw',X,Y),('rank_vehicle',X,Yvehicle),('value_raw',Xvalue,Y),('value_vehicle',Xvalue,Yvehicle)]
    if os.environ.get('YB_SPARSE_DIAGNOSTICS')=='1':
        variants=[('rank_raw',X,Y)]
        for k in [25,50,100,200]:
            selected=torch.zeros_like(X);selected.scatter_(1,torch.topk(X,k,dim=1).indices,1)
            selected-=selected.mean(1,keepdim=True);selected/=selected.std(1,keepdim=True,unbiased=False)
            variants.append((f'top{k}',selected,Y))
        (out/'sparse_plan.json').write_text(json.dumps({'question':'Does concentrating on the most relevant strains improve retrieval of known targets?',
             'candidate_sizes':[25,50,100,200],'response':'same frozen ranked barcode responses','calibration':'same 1000 shared-strain permutations and dose maximization',
             'selection':'only in training chemical-scaffold folds; best-known-target RR; dense baseline included','original_screen':'unchanged'},indent=2))
    for name,xx,yy in variants:
        observed=reduce_scores(xx@yy/S)
        sums=[torch.zeros((T,C),dtype=torch.float64,device=dev) for _ in range(2)]
        sums2=[torch.zeros_like(a) for a in sums]
        for seed in seeds:
            permutation_rng=np.random.default_rng(int(seed))
            perm=permutation_rng.permutation(S) if strata_labels is None else np.arange(S)
            if strata_labels is not None:
                for label in sorted(set(strata_labels)):
                    take=np.flatnonzero(strata_labels==label);perm[take]=permutation_rng.permutation(take)
            perm=torch.tensor(perm,device=dev)
            for k,a in enumerate(reduce_scores(xx[:,perm]@yy/S)):
                a=a.to(torch.float64);sums[k]+=a;sums2[k]+=a.square()
        for k,di in enumerate(['pos','neg']):
            mean=sums[k]/len(seeds);sd=torch.sqrt(torch.clamp(sums2[k]/len(seeds)-mean.square(),min=0))
            if not torch.all(sd>0):raise ValueError('degenerate diagnostic null')
            score=((observed[k].to(torch.float64)-mean)/sd).cpu().numpy()
            diagnostic[name+'_'+di]=score
        summary[name]={'task_transform':'centered standardized top-K binary membership' if name.startswith('top') else ('continuous ensemble score' if name.startswith('value') else 'rank'),
                       'response_transform':'vehicle median/MAD then rank' if name.endswith('vehicle') else 'rank of saved contrast',
                       'max_rho_cpu_difference':max(float(np.abs(a.cpu().numpy()-s[di+'_spearman_rho']).max()) for di,a in zip(['pos','neg'],observed)) if name=='rank_raw' else None}
        print('diagnostic completed',name,flush=True)
    np.savez_compressed(out/'diagnostic_scores.npz',**diagnostic)
    (out/'diagnostic_scoring_summary.json').write_text(json.dumps(summary,indent=2))
    sys.exit(0)

means={di:torch.tensor(s[di+'_null_mean'],dtype=torch.float32,device=dev) for di in ['pos','neg']}
sds={di:torch.tensor(s[di+'_null_sd'],dtype=torch.float32,device=dev) for di in ['pos','neg']}
obs_z={di:(s[di+'_spearman_rho']-s[di+'_null_mean'])/s[di+'_null_sd'] for di in ['pos','neg']}
obs_calc=reduce_scores(X@Y/S)
errors={di:float(np.max(np.abs(a.cpu().numpy()-s[di+'_spearman_rho']))) for di,a in zip(['pos','neg'],obs_calc)}
assert max(errors.values())<5e-7,errors  # float32 CPU/GPU accumulation tolerance
print('Archived score reproduction max errors:',errors,flush=True)

# Decomposition is performed after the actual screen's rank standardization.
mean_task=X.mean(0,keepdim=True)
common=(mean_task@Y/S).cpu().numpy()[0]
coverage_task=torch.zeros_like(mean_task)
if len(outside):
    coverage_task[:,support]=mean_task[:,support].mean()
    coverage_task[:,outside]=mean_task[:,outside].mean()
coverage=(coverage_task@Y/S).cpu().numpy()[0]
norm_total=float(X.square().sum().item())
decomp={'n_measured_strains':S,'n_represented_measured_strains':len(support),'n_fixed_tail_measured_strains':len(outside),
 'task_rank_norm_fraction_across_target_mean':float(T*mean_task.square().sum().item()/norm_total),
 'task_rank_norm_fraction_coverage_only':float(T*coverage_task.square().sum().item()/norm_total),
 'score_reproduction_max_absolute_error':errors,
 'interpretation':'Squared norm fractions of the rank-standardized screen input, distinct from pre-rank common-mode correction; not fractions of biological information.'}

nomination_record=json.loads((PANEL/'nomination_record.json').read_text())
allowlists=nomination_record.get('allowlists',{di:'pair_allowlist_'+di+'40.tsv' for di in ['pos','neg']})
nom=[]
source_to_local={int(c):i for i,c in enumerate(d['conditions'])}
for di in ['pos','neg']:
    f=PANEL/allowlists[di]
    frame=pd.read_csv(f,sep='\t')
    if nomination_record.get('selection_mode')=='fixed_reference':
        if os.environ.get('YB_SCREEN_MAX_ONLY')!='1':raise ValueError('fixed-reference comparison requires full-screen maximum tests')
        import hashlib
        reference=Path(nomination_record['reference_directory'])/allowlists[di]
        if hashlib.sha256(reference.read_bytes()).hexdigest()!=nomination_record['reference_hashes'][di] or reference.read_bytes()!=f.read_bytes():
            raise ValueError('fixed reference hypotheses changed')
        selected=([tp[t] for t in frame.target_id],[cp[k] for k in frame.inchikey])
    elif nomination_record.get('selection_mode')=='discovery_p':
        selected=np.where(s[di+'_emp_p']<=nomination_record['discovery_p_max'])
    else:
        selected=np.unravel_index(np.argsort(obs_z[di],axis=None)[-len(frame):],(T,C))
    archived_top=set((targets[i],compounds[j]) for i,j in zip(*selected))
    assert archived_top==set(zip(frame.target_id,frame.inchikey))
    for r in frame.itertuples():
        i,j=tp[r.target_id],cp[r.inchikey]
        k=source_to_local[int(s[di+'_source_condition_index'][i,j])]
        sg=1 if di=='pos' else -1
        nom.append({'direction':di,'target_id':r.target_id,'inchikey':r.inchikey,'ti':i,'ci':j,'rho':float(s[di+'_spearman_rho'][i,j]),'z':float(obs_z[di][i,j]),'coverage_component':float(sg*coverage[k]),'across_target_mean_component':float(sg*common[k]),'target_deviation':float(s[di+'_spearman_rho'][i,j]-sg*common[k])})
nom=pd.DataFrame(nom,columns=['direction','target_id','inchikey','ti','ci','rho','z','coverage_component','across_target_mean_component','target_deviation'])
nom.to_csv(OUT/'candidate_score_decomposition.tsv',sep='\t',index=False)
for di in ['pos','neg']:
    z=obs_z[di];sel=nom[nom.direction==di]
    idx=(sel.ti*C+sel.ci).to_numpy(dtype=int);cnt=np.bincount(idx%C,minlength=C)
    decomp['observed_'+di]={'n_nominees':len(sel),'nomination_min_z':float(sel.z.min()) if len(sel) else None,'nomination_max_z':float(sel.z.max()) if len(sel) else None,'unique_compounds':int((cnt>0).sum()),'largest_compound_count':int(cnt.max()),'mean_coverage_fraction_at_nominations':float((sel.coverage_component/sel.rho).mean()) if len(sel) else None,'mean_target_deviation':float(sel.target_deviation.mean()) if len(sel) else None}
(OUT/'decomposition_summary.json').write_text(json.dumps(decomp,indent=2))

if os.environ.get('YB_SCREEN_MAX_ONLY')=='1':
    output=OUT/'global_tail';output.mkdir(exist_ok=True)
    N=int(os.environ.get('YB_SCREEN_MAX_PERMUTATIONS','10000'))
    seed=int(os.environ.get('YB_SCREEN_MAX_SEED','20261015'))
    if N<1:raise ValueError('screen-wide permutations must be positive')
    rng=np.random.default_rng(seed);rows=[]
    (output/'analysis_plan.json').write_text(json.dumps(dict(n=N,seed=seed,n_targets=T,n_compounds=C,
        scope='Maximum over every target-compound pair and both directions; dose maximization repeated in every permutation',
        nomination_rule=nomination_record['rule'],n_directional_nominees=len(nom)),indent=2))
    for b in range(N):
        perm=torch.tensor(permute_strains(np.arange(S)),device=dev)
        maxima=[float(((a-means[di])/sds[di]).max()) for di,a in zip(['pos','neg'],reduce_scores(X[:,perm]@Y/S))]
        rows.append(dict(permutation=b+1,pos=maxima[0],neg=maxima[1],both=max(maxima)))
        if (b+1)%1000==0:print('full-screen maximum',b+1,flush=True)
    null=pd.DataFrame(rows);null.to_csv(output/'screen_max_null.tsv',sep='\t',index=False)
    pairs=[]
    for r in nom.itertuples():
        pairs.append(dict(direction=r.direction,target_id=r.target_id,inchikey=r.inchikey,observed_z=r.z,
            p_screen_max_same_direction=float((1+null[r.direction].ge(r.z).sum())/(N+1)),
            p_screen_max_both_directions=float((1+null.both.ge(r.z).sum())/(N+1))))
    pd.DataFrame(pairs,columns=['direction','target_id','inchikey','observed_z','p_screen_max_same_direction','p_screen_max_both_directions']).to_csv(output/'candidate_maxstat_tests.tsv',sep='\t',index=False)
    print('Screen-wide maximum tests complete for',len(pairs),'directional nominees.',flush=True);sys.exit(0)
if nomination_record.get('selection_mode')=='discovery_p':
    raise ValueError('Variable-count nominations require YB_SCREEN_MAX_ONLY; legacy top40 replay is inapplicable')

# Descriptive all-target context: selected target ranks are not validation.
cids={'SQFNCEKTEKNPCN-UHFFFAOYSA-N':'C1','CLQPEZWMMZFTBG-UHFFFAOYSA-N':'C2'}
context=[]
scope=pd.read_csv(UP_SCOPE:=Path('/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901/product_pipeline/product/transfer_route_b/inputs/universe_targets.tsv'),sep='\t')
for key,label in cids.items():
    j=cp[key];di='pos' if label=='C1' else 'neg';sg=1 if di=='pos' else -1
    nominated=set(nom.loc[(nom.inchikey==key)&(nom.direction==di),'target_id'])
    for i,t in enumerate(targets):
        k=source_to_local[int(s[di+'_source_condition_index'][i,j])]
        context.append({'compound':label,'inchikey':key,'target_id':t,'direction':di,'nominated':t in nominated,'rho':s[di+'_spearman_rho'][i,j],'z':obs_z[di][i,j],'coverage_component':sg*coverage[k],'across_target_mean_component':sg*common[k],'target_deviation':s[di+'_spearman_rho'][i,j]-sg*common[k]})
pd.DataFrame(context).to_csv(OUT/'C1_C2_all_target_scores.tsv',sep='\t',index=False)

if os.environ.get('YB_GLOBAL_TAIL_ONLY') == '1':
    # Re-select in every null screen: tests the aggregate nomination set,
    # not the significance of any individual selected pair.
    from product_execute_hiphop import bh_fdr
    output=OUT/'global_tail';output.mkdir(exist_ok=True)
    N=10000;seed=20261015;rng=np.random.default_rng(seed)
    metrics=['top40_min_z','top40_mean_z','unique_compounds']
    plan={'n':N,'seed':seed,'selection':'top40 re-selected per direction in each shared-strain null',
          'metrics':metrics,'directions':['pos','neg'],'multiplicity':'BH across all six direction-metric comparisons',
          'scope':'Exploratory aggregate screen-level evidence; does not replace individual maximum-statistic P or residual tests.'}
    (output/'analysis_plan.json').write_text(json.dumps(plan,indent=2))
    observed={}
    for di in ['pos','neg']:
        g=nom[nom.direction==di]
        observed[di]=[g.z.min(),g.z.mean(),g.inchikey.nunique()]
    rows=[]
    for b in range(N):
        perm=torch.tensor(permute_strains(np.arange(S)),device=dev)
        for di,a in zip(['pos','neg'],reduce_scores(X[:,perm]@Y/S)):
            z=(a-means[di])/sds[di];vals,idx=torch.topk(z.flatten(),40)
            rows.append({'permutation':b+1,'direction':di,'max_z':float(vals[0]),'top40_min_z':float(vals[-1]),
                         'top40_mean_z':float(vals.mean()),'unique_compounds':int(torch.unique(idx%C).numel())})
        if (b+1)%1000==0:print('global tail',b+1,flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(output/'selection_replay.tsv',sep='\t',index=False)
    results=[]
    for di in ['pos','neg']:
        null=frame[frame.direction==di]
        for k,m in enumerate(metrics):
            obs=float(observed[di][k]);v=null[m].to_numpy()
            results.append({'direction':di,'metric':m,'observed':obs,'null_mean':float(v.mean()),
                            'null_95_low':float(np.quantile(v,.025)),'null_95_high':float(np.quantile(v,.975)),
                            'p_higher':float((1+(v>=obs).sum())/(N+1))})
    result=pd.DataFrame(results);result['q_bh_six']=bh_fdr(result.p_higher.to_numpy())
    result.to_csv(output/'aggregate_tests.tsv',sep='\t',index=False)
    maxima=frame.pivot(index='permutation',columns='direction',values='max_z')
    jointmax=maxima.max(axis=1).to_numpy();pairtests=[]
    for r in nom.itertuples():
        pairtests.append(dict(direction=r.direction,target_id=r.target_id,inchikey=r.inchikey,observed_z=r.z,
            p_screen_max_same_direction=float((1+(maxima[r.direction].to_numpy()>=r.z).sum())/(N+1)),
            p_screen_max_both_directions=float((1+(jointmax>=r.z).sum())/(N+1))))
    pd.DataFrame(pairtests).to_csv(output/'candidate_maxstat_tests.tsv',sep='\t',index=False)
    print(result.to_string(index=False),flush=True);sys.exit(0)

N=1000
rng=np.random.default_rng(20261011)
for mode in ['represented_only','all_strains']:
    records=[];counts={di:np.zeros((N,C),dtype=np.int16) for di in ['pos','neg']}
    maxzs={di:np.empty(N) for di in ['pos','neg']}
    exceed={di:torch.zeros((T,C),device=dev,dtype=torch.int32) for di in ['pos','neg']}
    sums={di:torch.zeros((T,C),device=dev,dtype=torch.float64) for di in ['pos','neg']}
    sums2={di:torch.zeros((T,C),device=dev,dtype=torch.float64) for di in ['pos','neg']}
    observed={di:torch.tensor(s[di+'_spearman_rho'],device=dev,dtype=torch.float32) for di in ['pos','neg']}
    t0=time.time()
    for b in range(N):
        if mode=='represented_only':
            perm=torch.tensor(permute_strains(support),device=dev)
            q=Xsub[:,perm]@Ysub/S+fixed
        else:
            perm=torch.tensor(permute_strains(np.arange(S)),device=dev)
            q=X[:,perm]@Y/S
        for di,a in zip(['pos','neg'],reduce_scores(q)):
            exceed[di]+=(a>=observed[di]).to(torch.int32)
            sums[di]+=a.to(torch.float64);sums2[di]+=a.to(torch.float64).square()
            z=(a-means[di])/sds[di]
            vals,idx=torch.topk(z.flatten(),40)
            vals=vals.cpu().numpy();idx=idx.cpu().numpy();cnt=np.bincount(idx%C,minlength=C)
            counts[di][b]=cnt;maxzs[di][b]=float(vals[0])
            records.append({'permutation':b+1,'direction':di,'max_z':float(vals[0]),'top40_min_z':float(vals[-1]),'top40_mean_z':float(vals.mean()),'unique_compounds':int((cnt>0).sum()),'largest_compound_count':int(cnt.max()),'C1_count':int(cnt[cp[list(cids)[0]]]),'C2_count':int(cnt[cp[list(cids)[1]]])})
        if (b+1)%50==0:print(mode,b+1,'elapsed_s',round(time.time()-t0,1),flush=True)
    arrays={}
    for di in ['pos','neg']:
        arrays[di+'_null_mean']=(sums[di]/N).cpu().numpy()
        arrays[di+'_null_sd']=torch.sqrt(torch.clamp(sums2[di]/N-(sums[di]/N).square(),min=0)).cpu().numpy()
        arrays[di+'_emp_p']=((exceed[di].to(torch.float64)+1)/(N+1)).cpu().numpy()
        arrays[di+'_top40_counts']=counts[di]
    np.savez_compressed(OUT/(mode+'_null.npz'),**arrays)
    pd.DataFrame(records).to_csv(OUT/(mode+'_selection_replay.tsv'),sep='\t',index=False)
    result=[]
    jointmax=np.maximum(maxzs['pos'],maxzs['neg'])
    for r in nom.itertuples():
        result.append({'direction':r.direction,'target_id':r.target_id,'inchikey':r.inchikey,'observed_z':r.z,'p_fixed_pair':float(arrays[r.direction+'_emp_p'][r.ti,r.ci]),'p_screen_max_same_direction':float((1+(maxzs[r.direction]>=r.z).sum())/(N+1)),'p_screen_max_both_directions':float((1+(jointmax>=r.z).sum())/(N+1))})
    pd.DataFrame(result).to_csv(OUT/(mode+'_candidate_pvalues.tsv'),sep='\t',index=False)
    print('Finished',mode,flush=True)

(OUT/'mapping_controls_complete.json').write_text(json.dumps({'n':N,'seed':20261011,'modes':['represented_only','all_strains'],'device':str(torch.cuda.get_device_name(0)),'tf32':False,'analysis_dir':str(OUT),'nomination_dir':str(PANEL),'scope':'Retrospective sensitivity analyses on the supplied screens and nominations.'},indent=2))
