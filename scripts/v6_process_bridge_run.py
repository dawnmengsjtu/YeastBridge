#!/usr/bin/env python3
"""v6_process_bridge_run — V6 过程桥正式运行（B 骨架，降级 A 升级 B）。

协议：docs/cross_species_match/V6_PROCESS_BRIDGE_PROTOCOL.md
配置：configs/v6_process_bridge.json（输入 SHA-256 冻结，运行前强制核对）

部分：
  A  Norman 19 类升级融合（V5.1 轨道 B 同数据；成员 5→7；嵌套 (family,k)；
     jrepr=A 的 W_A 降级为可选成员）
  B  universe 靶点模块接力（76 靶点；LOO 签名；H/J/F 三臂 vs 置换零分布；门 R）
  C  应用层重接：B 模式任务表导出（导出分按预声明规则取 F 或 H）

用法：
  python scripts/v6_process_bridge_run.py --config configs/v6_process_bridge.json \
      --output <run_dir> --part all
"""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cross_species_match import load_matrix  # noqa: E402
from v5_joint_formal_run import (  # noqa: E402
    auc_p, build_projection_maps, l2n, project_v2_matrix, remove_pcs,
    sha256_file, signflip_p, z_rows)


# ---------------------------------------------------------------- setup
def v6_setup(config_path):
    cfg = json.loads(Path(config_path).read_text())
    root = Path(cfg["stage_root"])

    # 输入哈希核对（不匹配即中止）
    manifest = {}
    for k, v in cfg["inputs"].items():
        if not isinstance(v, dict) or "path" not in v:
            continue
        p = root / v["path"]
        h = sha256_file(str(p))
        manifest[k] = {"path": v["path"], "sha256": h, "match": h == v["sha256"]}
        if h != v["sha256"]:
            print(f"ABORT: input hash mismatch: {k} ({v['path']})")
            sys.exit(2)
    print("input hashes: all match", flush=True)

    import anndata as ad
    gate = cfg["data_quality_gate"]["column_abs_max_threshold"]

    # ---- K562（V5.1 main() 逐字口径） ----
    a = ad.read_h5ad(root / cfg["inputs"]["k562_h5ad"]["path"])
    genes_k = list(a.var.gene_name)
    obs_names = list(a.obs.index)
    pat = re.compile(r"^\d+_([^_]+)_(P1|P1P2)_")
    keep, target_of = [], []
    for name in obs_names:
        m = pat.match(name)
        if m and m.group(2) in {"P1", "P2", "P1P2"}:
            keep.append(name)
            target_of.append(m.group(1))
    Xp = a.X[[obs_names.index(n) for n in keep]]
    Xp = Xp.toarray() if hasattr(Xp, "toarray") else np.asarray(Xp)
    Xp = np.nan_to_num(Xp, nan=0.0)
    col_ok = np.abs(Xp).max(0) < gate
    n_drop = int((~col_ok).sum())
    assert n_drop == cfg["data_quality_gate"]["expected_dropped_columns"], n_drop
    Xp = Xp[:, col_ok]
    genes_k = [g for g, ok in zip(genes_k, col_ok) if ok]
    gidx_k = {g: j for j, g in enumerate(genes_k)}
    t2row = {}
    for j, t in enumerate(target_of):
        t2row.setdefault(t, j)
    print(f"K562 P1 profiles {Xp.shape[0]}, gate dropped {n_drop} columns", flush=True)
    del a

    # ---- 酵母矩阵 ----
    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    gidx_y = {g: j for j, g in enumerate(genes_y)}
    meta_y = pd.read_csv(root / cfg["inputs"]["yeast_metadata"]["path"], sep="\t")
    meta_y = meta_y[meta_y.benchmark_eligible.astype(str).str.lower() == "1"]
    meta_y = meta_y[meta_y.signature_id.isin(mat_y)].reset_index(drop=True)
    Y = np.stack([mat_y[s] for s in meta_y.signature_id])
    orf2row = {}
    for i, r in meta_y.iterrows():
        orf2row.setdefault(str(r.target_stable_id), i)
    del mat_y

    # ---- goSlim 投影（A 部分成员 + valid_k 口径） ----
    edges_gs, shared_gs = build_projection_maps(
        root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)
    V_k = project_v2_matrix(Xp, genes_k, "human", edges_gs, N_MOD)
    valid_k = np.linalg.norm(V_k, axis=1) > 1e-9
    V_y_gs = project_v2_matrix(Y, genes_y, "yeast", edges_gs, N_MOD)
    del V_k

    # ---- GO-BP ES 空间（V5.1 逐字口径 + 保留 B 稀疏阵） ----
    anc = pd.read_csv(root / cfg["inputs"]["go_bp_ancestors"]["path"],
                      sep="\t", compression="gzip",
                      usecols=["species", "gene", "ancestor_go_id"])
    tmin = cfg["representation"]["es_term_min_genes"]
    tmax = cfg["representation"]["es_term_max_genes"]

    def term_sets(df):
        d = {}
        for g, t in zip(df.gene, df.ancestor_go_id):
            d.setdefault(t, []).append(g)
        return {t: gs for t, gs in d.items() if tmin <= len(gs) <= tmax}

    ts_h = term_sets(anc[(anc.species == "human") & anc.gene.isin(gidx_k)])
    ts_y = term_sets(anc[(anc.species == "yeast") & anc.gene.isin(gidx_y)])
    del anc
    terms = sorted(set(ts_h) & set(ts_y))
    assert len(terms) == cfg["representation"]["expected_shared_terms"], len(terms)
    print(f"GO-BP shared terms: {len(terms)}", flush=True)

    def es_with_B(ts, idx_map, R):
        rows_, cols_ = [], []
        for j, t in enumerate(terms):
            for g in ts[t]:
                rows_.append(idx_map[g])
                cols_.append(j)
        B = sp.csr_matrix((np.ones(len(rows_), dtype=np.float32), (rows_, cols_)),
                          shape=(R.shape[1], len(terms)))
        root_arr = np.sqrt(np.array([len(ts[t]) for t in terms], dtype=np.float32))
        return l2n((z_rows(R) @ B) / root_arr[None, :]), root_arr, B

    E_k, root_h, B_h = es_with_B({t: ts_h[t] for t in terms}, gidx_k, Xp)
    E_y, _, _ = es_with_B({t: ts_y[t] for t in terms}, gidx_y, Y)
    E_k = E_k.astype(np.float64)
    E_y = E_y.astype(np.float64)
    pc = cfg["representation"]["pc_removal"]
    Ek1 = remove_pcs(E_k, pc)
    Ey1 = remove_pcs(E_y, pc)
    # E_k 的均值与 PC1 方向（LOO 单行去除用，与 remove_pcs 全矩阵口径一致）
    _Ec = E_k - E_k.mean(axis=0, keepdims=True)
    _, _, _Vt = np.linalg.svd(_Ec, full_matrices=False)
    ek_mean = E_k.mean(axis=0)
    ek_v0 = _Vt[:pc].copy()
    del _Ec, _Vt

    # ---- 邻域 ----
    nbr_df = pd.read_csv(root / cfg["inputs"]["neighborhood_terms"]["path"], sep="\t")
    h = sha256_file(str(root / cfg["inputs"]["neighborhood_terms"]["path"]))
    assert h == cfg["inputs"]["neighborhood_terms"]["sha256"], "neighborhood hash"
    nbr_terms = list(nbr_df.term_id)
    assert set(nbr_terms) <= set(terms), "neighborhood not subset of shared terms"
    nbr_idx = np.array([terms.index(t) for t in nbr_terms])
    assert len(nbr_idx) == cfg["neighborhood"]["n_terms"], len(nbr_idx)
    print(f"neighborhood terms: {len(nbr_idx)}", flush=True)

    # ---- 同源对（V5.1 main() 口径：valid_k 过滤 → 753） ----
    edges_df = pd.read_csv(root / cfg["inputs"]["orthodb_edges"]["path"], sep="\t")
    h_partner = {}
    for hh, yy in zip(edges_df.human_target_symbol, edges_df.yeast_deletion_target):
        h_partner.setdefault(hh, set()).add(yy)
    pairs = []
    for _, e in edges_df.iterrows():
        j, k2 = t2row.get(e.human_target_symbol), orf2row.get(e.yeast_deletion_target)
        if j is None or k2 is None or not valid_k[j]:
            continue
        pairs.append((e.human_target_symbol, e.yeast_deletion_target, e.split, j, k2))
    assert len(pairs) == 753, f"usable pairs {len(pairs)}, expected 753"
    ph = np.array([p[3] for p in pairs], dtype=int)
    py_ = np.array([p[4] for p in pairs], dtype=int)
    pair_tr = np.array([p[2] == "train" for p in pairs])

    # ---- W_A + PLS4（V5.1 armC 配方逐字；A 血统成员，V6 中降级为可选） ----
    ta_lambda = [1.0, 10.0, 100.0, 1000.0]
    ny = E_y.shape[0]
    rng = np.random.default_rng(20260911)  # V5.1 seed_row_negatives
    sym_of_row = {}
    for p in pairs:
        sym_of_row.setdefault(p[3], p[0])
    row_negs = {}
    for r in np.unique(ph):
        blk = {k2 for k2 in (orf2row.get(y) for y in h_partner.get(sym_of_row[r], ()))
               if k2 is not None}
        cand = np.setdiff1d(np.arange(ny), sorted(blk) or [-1])
        row_negs[r] = rng.choice(cand, size=min(30, len(cand)), replace=False)

    Xtr, Ytr = Ek1[ph[pair_tr]], Ey1[py_[pair_tr]]
    Mxy = Xtr.T @ Ytr
    best = (-1.0, None)
    for lam in ta_lambda:
        W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
        S = (Ek1 @ W) @ Ey1.T
        a_tr, _ = auc_p(S[ph[pair_tr], py_[pair_tr]],
                        np.concatenate([S[r, row_negs[r]] for r in np.unique(ph[pair_tr])]))
        if a_tr > best[0]:
            best = (a_tr, W)
    W_A, w_train_auc = best[1], best[0]
    mx_ = Xtr.mean(0); sx_ = np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_ = Ytr.mean(0); sy_ = np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                full_matrices=False)
    Wx4 = Usv[:, :4]; Wy4 = Vsv.T[:, :4]
    print(f"W_A refit (V5.1 recipe): ridge_train_auc={w_train_auc:.4f}", flush=True)

    # ---- Norman（V5.1 轨道 B 逐字口径） ----
    genes_n, mat_n = load_matrix(root / cfg["inputs"]["higher_matrix"]["path"])
    meta_n = pd.read_csv(root / cfg["inputs"]["higher_metadata"]["path"], sep="\t")
    meta_n = meta_n[meta_n.benchmark_eligible.astype(str).str.lower() == "true"]
    meta_n = meta_n[meta_n.signature_id.isin(mat_n)].reset_index(drop=True)
    assert len(meta_n) == 207, len(meta_n)
    Nn = np.stack([mat_n[s] for s in meta_n.signature_id])
    gidx_n = {g: j for j, g in enumerate(genes_n)}
    rows_n, cols_n = [], []
    for j_, t in enumerate(terms):
        for g in ts_h[t]:
            if g in gidx_n:
                rows_n.append(gidx_n[g])
                cols_n.append(j_)
    B_n = sp.csr_matrix((np.ones(len(rows_n), dtype=np.float32),
                         (rows_n, cols_n)), shape=(Nn.shape[1], len(terms)))
    E_n = l2n((z_rows(Nn) @ B_n) / root_h[None, :]).astype(np.float64)
    En1 = remove_pcs(E_n, pc)
    V_n = project_v2_matrix(Nn, genes_n, "human", edges_gs, N_MOD)
    del Nn, mat_n

    CLASSES = cfg["v6a_norman"].get("label_universe") or json.loads(
        (HERE.parent / "configs/v5_1_joint_formal.json").read_text())["trackB"]["label_universe"]
    CLS_IDX = {c: i for i, c in enumerate(CLASSES)}
    lab_n = meta_n.function_label.values
    grp_n = meta_n.split_group.values
    in_cls = np.array([str(l) in CLS_IDX for l in lab_n])
    yeast_ref_mask = np.array([(str(m.split).lower() == "train" and
                                str(m.function_label) in CLS_IDX)
                               for _, m in meta_y.iterrows()])

    return {
        "cfg": cfg, "root": root, "manifest": manifest,
        "Xp": Xp, "genes_k": genes_k, "gidx_k": gidx_k, "t2row": t2row,
        "meta_y": meta_y, "orf2row": orf2row, "Y": Y, "genes_y": genes_y,
        "terms": terms, "nbr_idx": nbr_idx, "nbr_terms": nbr_terms,
        "B_h": B_h, "root_h": root_h, "ek_mean": ek_mean, "ek_v0": ek_v0,
        "E_k": E_k, "E_y": E_y, "Ek1": Ek1, "Ey1": Ey1,
        "V_y_gs": V_y_gs, "V_n": V_n, "N_MOD": N_MOD,
        "W_A": W_A, "Wx4": Wx4, "Wy4": Wy4,
        "mx": mx_, "sx": sx_, "my": my_, "sy": sy_,
        "pairs": pairs, "ph": ph, "py": py_, "pair_tr": pair_tr,
        "h_partner": h_partner, "edges_df": edges_df,
        "meta_n": meta_n, "E_n": E_n, "En1": En1,
        "CLASSES": CLASSES, "CLS_IDX": CLS_IDX, "lab_n": lab_n,
        "grp_n": grp_n, "in_cls": in_cls, "yeast_ref_mask": yeast_ref_mask,
        "w_train_auc": w_train_auc,
    }


def cls_scores(qv, refs, ref_lab, n_cls):
    sims = refs @ qv
    out = np.full(n_cls, -np.inf)
    for c in range(n_cls):
        m = ref_lab == c
        if m.any():
            out[c] = sims[m].max()
    return out


def rank_of_truth(scores, truth, CLS_IDX):
    if truth not in CLS_IDX:
        return 0
    order = np.argsort(-scores, kind="stable")
    return int(np.where(order == CLS_IDX[truth])[0][0]) + 1


def label_centroids(mat, lab_idx, n_cls):
    return np.stack([mat[lab_idx == c].mean(0) if (lab_idx == c).any()
                     else np.zeros(mat.shape[1]) for c in range(n_cls)])


# ---------------------------------------------------------------- part A
def part_a(st, out_dir):
    cfgA = st["cfg"]["v6a_norman"]
    CLASSES, CLS_IDX = st["CLASSES"], st["CLS_IDX"]
    nC = len(CLASSES)
    meta_n, lab_n, grp_n, in_cls = st["meta_n"], st["lab_n"], st["grp_n"], st["in_cls"]
    En1, V_n = st["En1"], st["V_n"]
    Ey1, V_y_gs = st["Ey1"], st["V_y_gs"]
    W_A = st["W_A"]
    nbr_idx = st["nbr_idx"]
    yrm = st["yeast_ref_mask"]
    Yg0 = V_y_gs[yrm]
    Ey64r = st["E_y"][yrm]
    Ey1r = Ey1[yrm]
    lab_y_ref = np.array([CLS_IDX[str(m.function_label)]
                          for _, m in st["meta_y"][yrm].iterrows()])
    En1_nbr = l2n(En1[:, nbr_idx])
    Ey1r_nbr = l2n(Ey1r[:, nbr_idx])
    n = len(meta_n)

    # V5.1 五成员（逐字口径）
    S_h = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = En1[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h[i] = cls_scores(En1[i], Href, lab_h, nC)
    S_y_gs0 = np.stack([cls_scores(V_n[i], Yg0, lab_y_ref, nC) for i in range(n)])
    S_jrepr = np.stack([cls_scores(En1[i] @ W_A, Ey1r, lab_y_ref, nC) for i in range(n)])
    S_y_es0 = np.stack([cls_scores(st["E_n"][i], Ey64r, lab_y_ref, nC) for i in range(n)])
    lam_lab = 100.0  # V5.1 label_centroid_ridge_lambda
    S_lab = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        trm = (grp_n != g) & in_cls
        h_lab = np.array([CLS_IDX[str(x)] for x in lab_n[trm]])
        Ch = label_centroids(V_n[trm], h_lab, nC)
        Cy = label_centroids(Yg0, lab_y_ref, nC)
        W = np.linalg.solve(Ch.T @ Ch + lam_lab * np.eye(Ch.shape[1]), Ch.T @ Cy)
        for i in np.flatnonzero(grp_n == g):
            S_lab[i] = cls_scores(V_n[i] @ W, Yg0, lab_y_ref, nC)
    # 新邻域成员
    S_h_nbr = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = En1_nbr[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h_nbr[i] = cls_scores(En1_nbr[i], Href, lab_h, nC)
    S_y_nbr = np.stack([cls_scores(En1_nbr[i], Ey1r_nbr, lab_y_ref, nC) for i in range(n)])

    MEM = {"h_es1": S_h, "y_gs0": S_y_gs0, "jrepr": S_jrepr, "lab_gs0": S_lab,
           "y_es0": S_y_es0, "h_nbr": S_h_nbr, "y_nbr": S_y_nbr}
    families = {k: list(v) for k, v in cfgA["families"].items()}
    k_grid = list(cfgA["rrf_k_grid"])

    def rrf(members, k):
        out = np.zeros((n, nC))
        for m in members:
            S = MEM[m]
            for i in range(n):
                out[i] += 1.0 / (k + np.argsort(np.argsort(-S[i])))
        return out

    # 预计算全部 (family,k) 融合分与逐查询 rr
    rr_cfg = {}
    for fam, members in families.items():
        for k in k_grid:
            Sf = rrf(members, k)
            rr_cfg[(fam, k)] = np.array([
                1.0 / r if (r := rank_of_truth(Sf[i], str(lab_n[i]), CLS_IDX)) > 0 else 0.0
                for i in range(n)])
    rr_single = {}
    for m in MEM:
        S = MEM[m]
        rr_single[m] = np.array([
            1.0 / r if (r := rank_of_truth(S[i], str(lab_n[i]), CLS_IDX)) > 0 else 0.0
            for i in range(n)])

    single_fams = [f for f in families if f.startswith("single:")]

    def nested(cand_fams):
        rr_nest = np.full(n, np.nan)
        picks = {}
        for g in pd.unique(grp_n):
            te = grp_n == g
            tr = np.flatnonzero((~te) & in_cls)
            # 组内选 (family,k): MRR 最高; 平局规则(预注册): 成员数少者优先, 再 k 小者优先
            scored = []
            for fam in cand_fams:
                for k in k_grid:
                    v = float(rr_cfg[(fam, k)][tr].mean())
                    scored.append(((-v, len(families[fam]), k), (fam, k)))
            scored.sort(key=lambda x: x[0])
            fam_pick, k_pick = scored[0][1]
            picks[f"{fam_pick}|k={k_pick}"] = picks.get(f"{fam_pick}|k={k_pick}", 0) + 1
            for i in np.flatnonzero(te):
                rr_nest[i] = rr_cfg[(fam_pick, k_pick)][i]
        return rr_nest, picks

    rr_fusion, picks_fusion = nested(list(families))
    rr_comp, picks_comp = nested(single_fams)
    d = rr_fusion - rr_comp
    p = signflip_p(d[~np.isnan(d)], cfgA["seed_signflip_base"], cfgA["n_signflip"])
    gate = bool(p < 0.05)

    resA = {
        "fusion_mrr": float(np.nanmean(rr_fusion)),
        "comparator_best_single_mrr": float(np.nanmean(rr_comp)),
        "delta": float(np.nanmean(d)),
        "signflip_p": p,
        "gate_pass": gate,
        "fusion_picks": picks_fusion,
        "comparator_picks": picks_comp,
        "single_arm_mrr": {m: float(v.mean()) for m, v in rr_single.items()},
        "v51_reference": {"nested_joint_mrr": 0.348, "comparator_y_gs0": 0.286,
                          "signflip_p": 0.00244},
    }
    # 分层描述（强制功效注记）
    strat = {}
    for lab in cfgA["descriptive_strata"]:
        m = np.array([str(x) == lab for x in lab_n])
        if m.sum() == 0:
            continue
        strat[lab] = {"n": int(m.sum()),
                      "fusion_mrr": float(np.nanmean(rr_fusion[m])),
                      "comparator_mrr": float(np.nanmean(rr_comp[m])),
                      "delta": float(np.nanmean(d[m])),
                      "power_note": "descriptive only; n small — no gate claim"}
    resA["label_strata_descriptive"] = strat
    te_n = np.array([("test" in str(s).lower()) for s in meta_n.split])
    resA["official_test_descriptive"] = {
        "n": int(te_n.sum()),
        "fusion_mrr": float(np.nanmean(rr_fusion[te_n])),
        "comparator_mrr": float(np.nanmean(rr_comp[te_n]))}
    print(f"  A(V6) fusion MRR={resA['fusion_mrr']:.3f} comp={resA['comparator_best_single_mrr']:.3f} "
          f"Δ={resA['delta']:+.4f} p={p:.4g} gate={gate}", flush=True)

    pd.DataFrame({
        "signature_id": meta_n.signature_id, "label": lab_n, "split_group": grp_n,
        "official_test": te_n, "rr_comparator_best_single": rr_comp,
        "rr_nested_fusion": rr_fusion,
    }).to_csv(out_dir / "per_query_rr_trackA_v6.tsv", sep="\t", index=False)
    return resA


# ---------------------------------------------------------------- obo / GAF
def parse_obo_parents(obo_path):
    parents = {}
    cur_id, in_term = None, False
    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                in_term, cur_id = True, None
                continue
            if line.startswith("[") and line != "[Term]":
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cur_id = line[4:].strip()
                parents.setdefault(cur_id, set())
            elif cur_id is None:
                continue
            elif line.startswith("is_a: "):
                parents[cur_id].add(line[6:].split()[0].strip())
            elif line.startswith("relationship: part_of "):
                parents[cur_id].add(line[len("relationship: part_of "):].split()[0].strip())
    return parents


def ancestor_closure(term, parents, cache):
    if term in cache:
        return cache[term]
    seen, stack = {term}, [term]
    while stack:
        t = stack.pop()
        for p in parents.get(t, ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    cache[term] = seen
    return seen


def parse_gaf_bp(path, keys_by_col):
    """GAF 2.x 0-based 列: 2=DB_Object_Symbol, 3=Qualifier, 4=GO ID, 8=Aspect,
    10=DB_Object_Synonym(s)('|' 分隔, 酵母 ORF 在此列)。
    keys_by_col: {col_idx: set_of_keys}; 列 10 按 '|' 分词匹配。
    Returns {key: set(direct BP terms)}"""
    out = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or f[8] != "P" or "NOT" in f[3]:
                continue
            for col, keys in keys_by_col.items():
                if col >= len(f):
                    continue
                toks = f[col].split("|") if col == 10 else [f[col]]
                for k in toks:
                    if k in keys:
                        out.setdefault(k, set()).add(f[4])
    return out


# ---------------------------------------------------------------- part B
def part_b(st, out_dir):
    cfgB = st["cfg"]["v6b_relay"]
    root, cfg = st["root"], st["cfg"]
    nbr_terms, nbr_idx = st["nbr_terms"], st["nbr_idx"]
    nbr_pos = {t: i for i, t in enumerate(nbr_terms)}
    nN = len(nbr_terms)

    parents = parse_obo_parents(root / cfg["inputs"]["go_obo"]["path"])
    clo_cache = {}

    # ---- 靶点集合（冻结 C 臂导出记录 status=OK） ----
    rec = json.loads((root / cfg["inputs"]["c_arm_task_export_record"]["path"]).read_text())
    targets = [t for t, v in rec["targets"].items() if v.get("status") == "OK"]
    print(f"  B targets from frozen export record: {len(targets)}", flush=True)

    # ---- M_t：人 GAF（col2 = HGNC symbol） ----
    gaf_h = parse_gaf_bp(root / cfg["inputs"]["human_gaf"]["path"], {2: set(targets)})
    M = {}
    for t in targets:
        acc = set()
        for term in gaf_h.get(t, ()):
            acc |= (ancestor_closure(term, parents, clo_cache) & set(nbr_terms))
        M[t] = acc
    included = [t for t in targets if M[t]]
    dropped = [t for t in targets if not M[t]]
    print(f"  B M_t non-empty: {len(included)}, dropped: {len(dropped)} {dropped}", flush=True)
    if not included:
        raise SystemExit("FATAL: 无任何靶点获得非空 M_t —— 注释解析必然有误; "
                         "门 R 不予评估, 运行中止")

    # ---- a_s：酵母 GAF（ORF 本身 + 常用名双通道匹配, 同名列分词） ----
    orfs = [str(x) for x in st["meta_y"].target_stable_id]
    key2orf = {}
    for o in orfs:
        key2orf.setdefault(o, o)
    for o, sig in zip(orfs, st["meta_y"].signature_id):
        nm = str(sig).split(":")[-1].strip()
        if nm:
            key2orf.setdefault(nm, o)
    gaf_y = parse_gaf_bp(root / cfg["inputs"]["yeast_gaf"]["path"],
                         {2: set(key2orf), 10: set(key2orf)})
    terms_by_orf = {}
    for k, ts in gaf_y.items():
        terms_by_orf.setdefault(key2orf[k], set()).update(ts)
    A = np.zeros((len(orfs), nN))
    for i, o in enumerate(orfs):
        acc = set()
        for term in terms_by_orf.get(o, ()):
            acc |= (ancestor_closure(term, parents, clo_cache) & set(nbr_terms))
        for t_ in acc:
            A[i, nbr_pos[t_]] = 1.0
    n_ann = int((A.sum(1) > 0).sum())
    print(f"  B strains with nbr annotation: {n_ann}/{len(orfs)}", flush=True)
    A_n = l2n(A)

    Mmat = np.zeros((len(included), nN))
    for r, t in enumerate(included):
        for t_ in M[t]:
            Mmat[r, nbr_pos[t_]] = 1.0
    Mmat_n = l2n(Mmat)

    # ---- LOO 查询表示 ----
    Xp, gidx_k = st["Xp"], st["gidx_k"]
    B_h, root_h = st["B_h"], st["root_h"]
    ek_mean, ek_v0 = st["ek_mean"], st["ek_v0"]
    H_full = np.zeros((len(included), len(st["terms"])))
    for r, t in enumerate(included):
        j = st["t2row"][t]
        x = Xp[j].astype(np.float64).copy()
        if t in gidx_k:  # LOO: 靶点自身列置为行均值（预注册规则）
            x[gidx_k[t]] = x.mean()
        mu, sd = x.mean(), x.std()
        z = np.clip((x - mu) / (sd if sd > 0 else 1.0), -50, 50)
        e = z @ B_h / root_h
        ne = np.linalg.norm(e)
        e = e / ne if ne > 0 else e
        v = e - ek_mean
        for u in ek_v0:
            v = v - (v @ u) * u
        nv = np.linalg.norm(v)
        H_full[r] = v / nv if nv > 0 else v
    H_nbr = l2n(H_full[:, nbr_idx])
    Y_nbr = l2n(st["Ey1"][:, nbr_idx])
    s_H = Y_nbr @ H_nbr.T  # (strains, targets)

    # ---- J：冻结 C 臂导出分 verbatim ----
    task_dir = root / cfg["inputs"]["c_arm_task_dir"]
    s_J = np.full((len(orfs), len(included)), np.nan)
    for c, t in enumerate(included):
        fp = task_dir / f"yeast_task_{t}.tsv"
        df = pd.read_csv(fp, sep="\t")
        rowmap = {str(g): float(s) for g, s in zip(df.yeast_gene, df.ensemble_score)}
        s_J[:, c] = [rowmap.get(o, np.nan) for o in orfs]
    assert not np.isnan(s_J).any(), "J scores incomplete vs frozen export"
    gated = np.array([len(st["h_partner"].get(t, ())) > 0 for t in included])
    print(f"  B J-arm gating (has orthodb partner): {int(gated.sum())}/{len(included)}", flush=True)

    # ---- F：RRF k=10（J 门控参与） ----
    def ranks_desc(col):
        order = np.argsort(-col, kind="stable")
        rk = np.empty(len(col), dtype=np.int64)
        rk[order] = np.arange(1, len(col) + 1)
        return rk

    kR = 10
    s_F = np.empty_like(s_H)
    for c in range(len(included)):
        rh = 1.0 / (kR + ranks_desc(s_H[:, c]))
        if gated[c]:
            rj = 1.0 / (kR + ranks_desc(s_J[:, c]))
            s_F[:, c] = rh + rj
        else:
            s_F[:, c] = rh

    # ---- 一致性 + 置换 ----
    topk = cfgB["top_k_strains"]
    def top_mean_concord(S):
        idx = np.argsort(-S, axis=0, kind="stable")[:topk]      # (topk, targets)
        M50 = A_n[idx].mean(axis=0)                              # (targets_prof, nN)
        return M50 @ Mmat_n.T                                   # (prof, gold_target)

    CM_H = top_mean_concord(s_H)
    CM_J = top_mean_concord(s_J)
    CM_F = top_mean_concord(s_F)
    n_t = len(included)
    obs = {arm: float(np.mean(np.diag(CM))) for arm, CM in
           (("H", CM_H), ("J", CM_J), ("F", CM_F))}
    # 次统计量：top50 命中非空占比
    def top_hit_rate(S):
        idx = np.argsort(-S, axis=0, kind="stable")[:topk]
        hits = (A[idx] @ Mmat.T) > 0        # (topk, prof, gold)? 内存: 50*76*76 小
        return hits.mean(axis=0)            # (prof, gold)
    HR_H, HR_J, HR_F = top_hit_rate(s_H), top_hit_rate(s_J), top_hit_rate(s_F)
    obs2 = {arm: float(np.mean(np.diag(HR))) for arm, HR in
            (("H", HR_H), ("J", HR_J), ("F", HR_F))}

    for arm_, v_ in obs.items():
        assert np.isfinite(v_), f"FATAL: obs[{arm_}] non-finite"
    n_perm = cfgB["n_perm"]
    rng = np.random.default_rng(cfgB["seed_perm"])
    perm_stats = {arm: np.empty(n_perm) for arm in ("H", "J", "F")}
    perm_stats2 = {arm: np.empty(n_perm) for arm in ("H", "J", "F")}
    for b in range(n_perm):
        pi = rng.permutation(n_t)
        for arm, CM, HR in (("H", CM_H, HR_H), ("J", CM_J, HR_J), ("F", CM_F, HR_F)):
            perm_stats[arm][b] = np.mean(CM[pi, np.arange(n_t)])
            perm_stats2[arm][b] = np.mean(HR[pi, np.arange(n_t)])
    pvals, pvals2 = {}, {}
    for arm in ("H", "J", "F"):
        pvals[arm] = float((1 + (perm_stats[arm] >= obs[arm]).sum()) / (1 + n_perm))
        pvals2[arm] = float((1 + (perm_stats2[arm] >= obs2[arm]).sum()) / (1 + n_perm))
    gate_R = bool(pvals["H"] < 0.05 or pvals["F"] < 0.05)

    resB = {
        "n_targets_included": n_t, "n_targets_dropped_empty_Mt": dropped,
        "n_strains_annotated": n_ann, "j_gated_in": int(gated.sum()),
        "primary_stat_mean_top50_cos_concordance": obs,
        "primary_perm_p": pvals,
        "secondary_stat_top50_hit_rate": obs2,
        "secondary_perm_p": pvals2,
        "gate_R_pass": gate_R,
        "gated_targets": [t for t, g in zip(included, gated) if g],
    }
    print(f"  B relay: obs={ {k: round(v,4) for k,v in obs.items()} } "
          f"p={ {k: round(v,4) for k,v in pvals.items()} } gateR={gate_R}", flush=True)

    pd.DataFrame({
        "target": included,
        "Mt_terms": ["|".join(sorted(M[t])) for t in included],
        "obs_H": np.diag(CM_H), "obs_J": np.diag(CM_J), "obs_F": np.diag(CM_F),
        "hit_H": np.diag(HR_H), "hit_J": np.diag(HR_J), "hit_F": np.diag(HR_F),
        "j_gated": gated,
    }).to_csv(out_dir / "per_target_concordance.tsv", sep="\t", index=False)
    np.save(out_dir / "scores_H.npy", s_H)
    np.save(out_dir / "scores_F.npy", s_F)
    # H/F 全分数供 C 部分复用（含全部 included 靶点）
    (out_dir / "relay_included_targets.json").write_text(json.dumps(included))
    return resB, {"included": included, "s_H": s_H, "s_F": s_F,
                  "gated": gated, "H_full": H_full, "gate_R": gate_R,
                  "A": A, "orfs": orfs}


# ---------------------------------------------------------------- part C
def part_c(st, out_dir, relay):
    cfgC = st["cfg"]["v6c_export"]
    root = st["root"]
    included = relay["included"]
    gate_R = relay["gate_R"]
    use = "F" if gate_R else "H"   # 预声明决策规则, 无其他分支
    S = relay["s_F"] if gate_R else relay["s_H"]
    A = relay["A"]
    orfs = relay["orfs"]
    nbr_mask = A.sum(1) > 0
    Ey1 = st["Ey1"]

    # 尾段排序分：全 639 维 cos(h_t_full, Ey1 行)
    H_full_n = l2n(relay["H_full"])
    tail_raw = Ey1 @ H_full_n.T  # (strains, targets)

    task_dir = root / cfgC["outdir"]
    task_dir.mkdir(parents=True, exist_ok=True)
    report = {"mode": "export_gobridge_tasks_v6",
              "protocol": "V6_PROCESS_BRIDGE_PROTOCOL.md",
              "config_sha256": sha256_file(str(root / "yeastbridge_re_mvp/configs/v6_process_bridge.json")),
              "script_sha256": sha256_file(str(Path(__file__).resolve())),
              "decision_rule": "gate_R pass -> F else H",
              "gate_R": gate_R, "export_arm": use,
              "pool": len(orfs), "nbr_strains": int(nbr_mask.sum()),
              "targets": {}}
    for c, t in enumerate(included):
        scores = np.full(len(orfs), -np.inf)
        nb = S[nbr_mask, c]
        scores[nbr_mask] = nb
        min_nbr = float(nb.min())
        tl_idx = np.flatnonzero(~nbr_mask)
        tl = tail_raw[tl_idx, c]
        order_tl = tl_idx[np.argsort(-tl, kind="stable")]
        # 单调拼接：尾段严格低于邻域段最低分
        scores[order_tl] = min_nbr - 0.001 - np.arange(len(order_tl)) * 1e-9
        order = np.argsort(-scores, kind="stable")
        fp = task_dir / f"yeast_task_{t}.tsv"
        with fp.open("w", encoding="utf-8", newline="") as fh:
            fh.write("rank\tyeast_gene\tensemble_score\n")
            for rk, i in enumerate(order, 1):
                fh.write(f"{rk}\t{orfs[i]}\t{scores[i]:.9f}\n")
        report["targets"][t] = {"status": "OK", "rows": len(orfs), "path": str(fp)}
    (task_dir / "export_record.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"  C export arm={use} targets={len(included)} -> {task_dir}", flush=True)
    return report


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--part", default="all", choices=["a", "b", "c", "all"])
    args = ap.parse_args()

    out_dir = Path(args.output)
    (out_dir / "trackA_norman").mkdir(parents=True, exist_ok=True)
    (out_dir / "relay").mkdir(parents=True, exist_ok=True)

    st = v6_setup(args.config)
    res = {"protocol": st["cfg"]["protocol"],
           "config_sha256": sha256_file(args.config),
           "script_sha256": sha256_file(str(Path(__file__).resolve()))}
    (out_dir / "INPUT_MANIFEST.json").write_text(json.dumps(st["manifest"], indent=1))

    if args.part in ("a", "all"):
        res["v6a_norman"] = part_a(st, out_dir / "trackA_norman")
        (out_dir / "result_partial.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))

    relay_state = None
    if args.part in ("b", "all"):
        resB, relay_state = part_b(st, out_dir / "relay")
        res["v6b_relay"] = resB
        (out_dir / "result_partial.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))

    if args.part in ("c", "all"):
        if relay_state is None:
            # 从部分结果恢复（c 单独跑时）
            inc = json.loads((out_dir / "relay/relay_included_targets.json").read_text())
            s_H = np.load(out_dir / "relay/scores_H.npy")
            s_F = np.load(out_dir / "relay/scores_F.npy")
            gate_R = res.get("v6b_relay", {}).get("gate_R_pass")
            if gate_R is None:
                rp = json.loads((out_dir / "result_partial.json").read_text())
                gate_R = rp["v6b_relay"]["gate_R_pass"]
            # H_full 与 A 需要重建：重跑 setup 的 LOO/GAF 太重, c 单独跑不支持——要求 b+c 连跑
            raise SystemExit("part c 需要与 b 连跑（--part all 或 --part b 后同进程 c）")
        rep = part_c(st, out_dir, relay_state)
        res["v6c_export"] = {"export_arm": rep["export_arm"], "gate_R": rep["gate_R"],
                             "outdir": rep["outdir"] if "outdir" in rep else str(Path(st["cfg"]["v6c_export"]["outdir"]))}

    (out_dir / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("saved:", out_dir / "result.json")
    a = res.get("v6a_norman", {})
    b = res.get("v6b_relay", {})
    print(f"V6 gate B'(norman fusion): {a.get('gate_pass')} | gate R(relay): {b.get('gate_R_pass')}")


if __name__ == "__main__":
    main()
