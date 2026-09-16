#!/usr/bin/env python3
"""v7_b2_formal_run — V7 B2 轴联合增益正式运行（三臂两门，V5.1 同机换表示）。

协议：docs/cross_species_match/V7_B2_JOINT_FORMAL_PROTOCOL.md
配置：configs/v7_b2_joint_formal.json（输入 SHA-256 冻结，运行前强制核对）

用法：
  python scripts/v7_b2_formal_run.py --config configs/v7_b2_joint_formal.json \
      --output <run_dir>
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cross_species_match import load_matrix  # noqa: E402
from v5_joint_formal_run import (  # noqa: E402
    auc_p, build_projection_maps, l2n, project_v2_matrix, remove_pcs,
    sha256_file, signflip_p)


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


def pct_rows(S):
    M = np.empty_like(S)
    for r in range(S.shape[0]):
        M[r] = (np.argsort(np.argsort(S[r])) + 1) / S.shape[1]
    return M


# ------------------------------------------------------------------ setup
def v7_setup(config_path):
    cfg = json.loads(Path(config_path).read_text())
    root = Path(cfg["stage_root"])
    manifest = {}
    for k, v in cfg["inputs"].items():
        if not isinstance(v, dict):
            continue
        for path_key, sha_key in (("path", "sha256"), ("index_path", "index_sha256"),
                                  ("path_abs", "sha256"), ("index_path_abs", "index_sha256"),
                                  ("gene_order_path", "gene_order_sha256")):
            if path_key not in v:
                continue
            p = Path(v[path_key])
            p = p if p.is_absolute() else root / p
            h = sha256_file(str(p))
            want = v.get(sha_key)
            manifest[f"{k}.{path_key}"] = {"path": str(p), "sha256": h,
                                           "match": (h == want) if want else None}
            if want and h != want:
                print(f"ABORT: input hash mismatch: {k}.{path_key} ({p})")
                sys.exit(2)
    print("input hashes: all match", flush=True)

    def rp(rel):
        p = Path(rel)
        return p if p.is_absolute() else root / rel

    # ---- 人侧 ESM2 ----
    E_h = np.load(rp(cfg["inputs"]["esm2_human"]["path"])).astype(np.float64)
    idx_h = pd.read_csv(rp(cfg["inputs"]["esm2_human"]["index_path"]), sep="\t").fillna("")
    sym2row_h = {}
    for i, s in enumerate(idx_h["common"].astype(str)):
        sym2row_h.setdefault(s, i)
    print(f"human esm2: {E_h.shape}, symbols {len(sym2row_h)}", flush=True)

    # ---- B2 注入投影（esm2_cross_match.route_b_top 逐字同式） ----
    import torch
    ck = torch.load(rp(cfg["inputs"]["route_b_model"]["path"]),
                    map_location="cpu", weights_only=False)
    W_inj = ck["state_dict"]["pos_emb.proj.weight"].float().numpy().astype(np.float64)
    b_inj = ck["state_dict"]["pos_emb.proj.bias"].float().numpy().astype(np.float64)
    del ck

    # ---- 酵母侧三表 ----
    T_rb = np.load(rp(cfg["inputs"]["route_b_table"]["path"])).astype(np.float64)
    order_rb = pd.read_csv(rp(cfg["inputs"]["route_b_gene_order"]["path"]),
                           sep="\t", dtype=str)["systematic"].tolist()
    rb_row = {g: i for i, g in enumerate(order_rb)}
    E_y = np.load(rp(cfg["inputs"]["yeast_esm2"]["path_abs"])).astype(np.float64)
    idx_y = pd.read_csv(rp(cfg["inputs"]["yeast_esm2"]["index_path_abs"]), sep="\t")
    ye_row = {g: i for i, g in enumerate(idx_y["systematic"].astype(str))}
    G_scy = np.load(rp(cfg["inputs"]["scyeast_prior"]["path"])).astype(np.float64)
    scy_order = [l.strip() for l in open(rp(cfg["inputs"]["scyeast_prior"]["gene_order_path"]))]
    scy_row = {g: i for i, g in enumerate(scy_order)}

    # ---- 池 ----
    meta_y = pd.read_csv(root / cfg["inputs"]["yeast_metadata"]["path"], sep="\t")
    meta_y = meta_y[meta_y.benchmark_eligible.astype(str).str.lower() == "1"].reset_index(drop=True)
    orfs_all = [str(x) for x in meta_y.target_stable_id]
    pool = [o for o in orfs_all if o in rb_row and o in ye_row]
    pool_pos = {o: i for i, o in enumerate(pool)}
    n_pool = len(pool)
    scy_pool = [o for o in pool if o in scy_row]
    scy_pos = {o: i for i, o in enumerate(scy_pool)}
    print(f"pool: {n_pool}/{len(orfs_all)} (scyeast-covered {len(scy_pool)})", flush=True)

    Y_rb = l2n(T_rb[[rb_row[o] for o in pool]])
    Y_es = l2n(E_y[[ye_row[o] for o in pool]])
    Y_scy = l2n(G_scy[[scy_row[o] for o in scy_pool]])
    X_es_all = l2n(E_h)
    X_b2_all = l2n(E_h @ W_inj.T + b_inj)
    X_es1_all = remove_pcs(X_es_all, 1)
    X_b21_all = remove_pcs(X_b2_all, 1)
    Y_rb1 = remove_pcs(Y_rb, 1)
    Y_es1 = remove_pcs(Y_es, 1)

    # ---- 对集 ----
    edges_df = pd.read_csv(root / cfg["inputs"]["orthodb_edges"]["path"], sep="\t")
    h_partner = {}
    for hh, yy in zip(edges_df.human_target_symbol, edges_df.yeast_deletion_target):
        h_partner.setdefault(hh, set()).add(yy)
    v51 = pd.read_csv(root / cfg["inputs"]["v51_pair_list"]["path"], sep="\t")
    assert len(v51) == 753, len(v51)
    prim = [(str(r.human), str(r.yeast_orf), str(r.split), int(r.fold))
            for _, r in v51.iterrows()
            if str(r.human) in sym2row_h and str(r.yeast_orf) in pool_pos]
    print(f"primary pairs (V5.1 frozen ∩ B2-available): {len(prim)}/753", flush=True)
    prim_set = {(h, y) for h, y, _, _ in prim}
    ext = [(str(e.human_target_symbol), str(e.yeast_deletion_target), str(e.split))
           for _, e in edges_df.iterrows()
           if (str(e.human_target_symbol), str(e.yeast_deletion_target)) not in prim_set
           and str(e.human_target_symbol) in sym2row_h
           and str(e.yeast_deletion_target) in pool_pos]
    print(f"extended pairs (descriptive): {len(ext)}", flush=True)

    # ---- row_negs（V5.1 同配方, 池内重抽, 迭代序=查询符号排序） ----
    query_genes = sorted({h for h, _, _, _ in prim})
    q_pos = {g: i for i, g in enumerate(query_genes)}
    seed_neg = cfg["trackA"]["seed_row_negatives"]
    row_negs = {}
    rng = np.random.default_rng(seed_neg)
    for g in query_genes:
        blk = sorted(pool_pos[y] for y in h_partner.get(g, ()) if y in pool_pos)
        cand = np.setdiff1d(np.arange(n_pool), blk or [-1])
        row_negs[g] = rng.choice(cand, size=min(30, len(cand)), replace=False)

    return {"cfg": cfg, "root": root, "manifest": manifest,
            "sym2row_h": sym2row_h, "W_inj": W_inj, "b_inj": b_inj,
            "X_es_all": X_es_all, "X_b2_all": X_b2_all,
            "X_es1_all": X_es1_all, "X_b21_all": X_b21_all,
            "pool": pool, "pool_pos": pool_pos, "n_pool": n_pool,
            "Y_rb": Y_rb, "Y_es": Y_es, "Y_rb1": Y_rb1, "Y_es1": Y_es1,
            "Y_scy": Y_scy, "scy_pool": scy_pool, "scy_pos": scy_pos,
            "prim": prim, "ext": ext, "query_genes": query_genes, "q_pos": q_pos,
            "row_negs": row_negs, "h_partner": h_partner,
            "meta_y": meta_y, "edges_df": edges_df}


# ------------------------------------------------------------------ track A
def fit_ridge_lambda(Xq, Y, q_rows, y_rows, row_negs, query_genes, lam_grid):
    """V5.1 配方: λ 网格, 按给定配对 train AUC 选参（pos=配对行, neg=唯一查询行
    的 30 株冻结负例）。返回 (train_auc, W)。"""
    Xtr, Ytr = Xq[q_rows], Y[y_rows]
    Mxy = Xtr.T @ Ytr
    uniq_q = np.unique(q_rows)
    best = (-1.0, None)
    for lam in lam_grid:
        W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
        S = (Xq @ W) @ Y.T
        a_tr, _ = auc_p(S[q_rows, y_rows],
                        np.concatenate([S[qq, row_negs[query_genes[qq]]]
                                        for qq in uniq_q]))
        if a_tr > best[0]:
            best = (a_tr, W)
    return best


def track_a_armset(st, Xall, X1all, Y, Y1, label, out):
    """同构三臂机器: A 无对齐 / B 去PC1 / C ridge+PLS4 秩集成。Xall 行=全部人蛋白,
    查询切片按 st['query_genes']。"""
    cfgA = st["cfg"]["trackA"]
    prim = st["prim"]
    q_pos, pool_pos = st["q_pos"], st["pool_pos"]
    row_negs, query_genes = st["row_negs"], st["query_genes"]
    n_pool = Y.shape[0]

    hrow = [st["sym2row_h"][g] for g in query_genes]
    Xq, Xq1 = Xall[hrow], X1all[hrow]
    q_rows = np.array([q_pos[h] for h, _, _, _ in prim])
    y_rows = np.array([pool_pos[y] for _, y, _, _ in prim])
    split_tr = np.array([s == "train" for _, _, s, _ in prim])
    split_te = np.array([s == "test" for _, _, s, _ in prim])
    folds = np.array([f for _, _, _, f in prim])

    S_A = Xq @ Y.T
    S_B = Xq1 @ Y1.T

    def pls4_scores(Xtr_, Ytr_):
        mx_, sx_ = Xtr_.mean(0), np.where(Xtr_.std(0) > 0, Xtr_.std(0), 1)
        my_, sy_ = Ytr_.mean(0), np.where(Ytr_.std(0) > 0, Ytr_.std(0), 1)
        U, _, Vt = np.linalg.svd(((Xtr_ - mx_) / sx_).T @ ((Ytr_ - my_) / sy_),
                                 full_matrices=False)
        Px = l2n(((Xq - mx_) / sx_) @ U[:, :4])
        Py = l2n(((Y - my_) / sy_) @ Vt[:4].T)
        return Px @ Py.T

    tr_auc, W = fit_ridge_lambda(Xq, Y, q_rows[split_tr], y_rows[split_tr],
                                 row_negs, query_genes, cfgA["ridge_lambda_grid"])
    S_ridge = (Xq @ W) @ Y.T
    S_pls = pls4_scores(Xq[q_rows[split_tr]], Y[y_rows[split_tr]])
    S_C = np.mean([pct_rows(S_ridge), pct_rows(S_pls)], axis=0)

    resA = {"repr": label, "ridge_train_auc": float(tr_auc),
            "n_pairs_primary": len(prim), "n_pool": n_pool}
    te_uq = np.unique(q_rows[split_te])

    def arm_eval(S, name):
        pos = S[q_rows[split_te], y_rows[split_te]]
        neg = np.concatenate([S[qq, row_negs[query_genes[qq]]] for qq in te_uq])
        a, p = auc_p(pos, neg)
        resA[f"{name}_test_auc"] = a
        resA[f"{name}_test_p"] = p
        print(f"  A[{label}] {name}: test auc={a:.4f} p={p:.2g}", flush=True)

    for nm, S in (("armA_no_align", S_A), ("armB_unsup_pc1", S_B), ("armC_joint", S_C)):
        arm_eval(S, nm)

    def boot_delta(Sa, Sb, seed):
        # V5.1 逐字: 池化 pos/neg, 逐对 win-rate, bootstrap 重采样对
        pos_a = Sa[q_rows[split_te], y_rows[split_te]]
        neg_a = np.concatenate([Sa[qq, row_negs[query_genes[qq]]] for qq in te_uq])
        pos_b = Sb[q_rows[split_te], y_rows[split_te]]
        neg_b = np.concatenate([Sb[qq, row_negs[query_genes[qq]]] for qq in te_uq])
        wa = (pos_a[:, None] > neg_a[None, :]).mean(axis=1)
        wb = (pos_b[:, None] > neg_b[None, :]).mean(axis=1)
        rb_ = np.random.default_rng(seed)
        idx = np.arange(len(wa))
        ds = np.empty(cfgA["n_bootstrap"])
        for t in range(cfgA["n_bootstrap"]):
            s = rb_.choice(idx, len(idx), replace=True)
            ds[t] = wa[s].mean() - wb[s].mean()
        return float(ds.mean()), float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))

    for nm, Sa, Sb in (("C_vs_A", S_C, S_A), ("C_vs_B", S_C, S_B)):
        m, lo, hi = boot_delta(Sa, Sb, cfgA["seed_bootstrap"])
        resA[f"boot_{nm}"] = {"mean": m, "ci95_low": lo, "ci95_high": hi}
        print(f"  A[{label}] ΔAUC {nm}: {m:+.4f} CI[{lo:+.4f},{hi:+.4f}]", flush=True)

    # OOF（fold 标签复用 V5.1 冻结清单）
    oof = {nm: np.full(len(prim), np.nan) for nm in ("A", "B", "C")}
    oof_hit = {nm: np.full(len(prim), np.nan) for nm in ("B", "C")}

    def pct_ranks(S, mask):
        return np.array([(np.sum(S[qq] > S[qq, yy]) + 1) / n_pool
                         for qq, yy in zip(q_rows[mask], y_rows[mask])])

    for f in sorted(set(folds.tolist())):
        m_tr, m_te = folds != f, folds == f
        _, Wf = fit_ridge_lambda(Xq, Y, q_rows[m_tr], y_rows[m_tr],
                                 row_negs, query_genes, cfgA["ridge_lambda_grid"])
        Sc1 = (Xq @ Wf) @ Y.T
        Sp = pls4_scores(Xq[q_rows[m_tr]], Y[y_rows[m_tr]])
        SC_f = np.mean([pct_rows(Sc1), pct_rows(Sp)], axis=0)
        for nm, S in (("A", S_A), ("B", S_B), ("C", SC_f)):
            oof[nm][m_te] = pct_ranks(S, m_te)
        for i in np.flatnonzero(m_te):
            g = query_genes[q_rows[i]]
            partners = [pool_pos[y2] for y2 in st["h_partner"].get(g, ()) if y2 in pool_pos]
            if not partners:
                continue
            oof_hit["B"][i] = min(int(np.sum(S_B[q_rows[i]] > S_B[q_rows[i], c_])) + 1
                                  for c_ in partners) <= 10
            oof_hit["C"][i] = min(int(np.sum(SC_f[q_rows[i]] > SC_f[q_rows[i], c_])) + 1
                                  for c_ in partners) <= 10
    valid = ~np.isnan(oof["A"])
    sf_seed = cfgA["seed_signflip_base"]
    for nm_a, nm_b in (("C", "A"), ("C", "B")):
        d = oof[nm_b][valid] - oof[nm_a][valid]
        sf_seed += 1
        p_ = signflip_p(d, sf_seed, cfgA["n_signflip"])
        resA[f"oof_{nm_a}_vs_{nm_b}"] = {"gain": float(d.mean()), "signflip_p": p_}
        print(f"  A[{label}] OOF {nm_a}-vs-{nm_b}: {d.mean():+.4f} p={p_:.3g}", flush=True)
    hb = oof_hit["B"][~np.isnan(oof_hit["B"])]
    hc = oof_hit["C"][~np.isnan(oof_hit["C"])]
    if len(hb) and len(hc):
        d_hit = hc.astype(float) - hb.astype(float)
        sf_seed += 1
        resA["any_hit10"] = {"rate_B": float(hb.mean()), "rate_C": float(hc.mean()),
                             "signflip_p": signflip_p(d_hit, sf_seed, cfgA["n_signflip"])}
    else:
        resA["any_hit10"] = {"rate_B": None, "rate_C": None, "signflip_p": 1.0}
    print(f"  A[{label}] any-hit@10: B={resA['any_hit10']['rate_B']} "
          f"C={resA['any_hit10']['rate_C']} p={resA['any_hit10']['signflip_p']:.3g}", flush=True)

    resA["gate_primary"] = bool(resA["armC_joint_test_p"] < 0.05 and
                                resA["boot_C_vs_A"]["ci95_low"] > 0)
    resA["gate_secondary"] = bool(resA["oof_C_vs_A"]["signflip_p"] < 0.05 or
                                  resA["any_hit10"]["signflip_p"] < 0.05)
    print(f"  A[{label}] gates: primary={resA['gate_primary']} "
          f"secondary={resA['gate_secondary']}", flush=True)

    pd.DataFrame({"human": [h for h, _, _, _ in prim],
                  "yeast_orf": [y for _, y, _, _ in prim],
                  "split": [s for _, _, s, _ in prim], "fold": folds,
                  "oof_pct_rank_A": oof["A"], "oof_pct_rank_B": oof["B"],
                  "oof_pct_rank_C": oof["C"]}).to_csv(
        out / f"per_query_ranks_trackA_{label}.tsv", sep="\t", index=False)
    return resA, W


def aligner_cell(st, Xall, Y, orf_order, label):
    """附加格: 仅 ridge 对齐器 (repr_grid 式), 负例在该格候选空间内同配方重抽。"""
    cfgA = st["cfg"]["trackA"]
    pos_of = {o: i for i, o in enumerate(orf_order)}
    prim = [(h, y, s, f) for h, y, s, f in st["prim"] if y in pos_of]
    q_rows = np.array([st["q_pos"][h] for h, _, _, _ in prim])
    y_rows = np.array([pos_of[y] for _, y, _, _ in prim])
    split_tr = np.array([s == "train" for _, _, s, _ in prim])
    split_te = np.array([s == "test" for _, _, s, _ in prim])
    hrow = [st["sym2row_h"][g] for g in st["query_genes"]]
    Xq = Xall[hrow]

    rng = np.random.default_rng(cfgA["seed_row_negatives"])
    row_negs = {}
    for g in st["query_genes"]:
        blk = sorted(pos_of[y] for y in st["h_partner"].get(g, ()) if y in pos_of)
        cand = np.setdiff1d(np.arange(len(orf_order)), blk or [-1])
        row_negs[g] = rng.choice(cand, size=min(30, len(cand)), replace=False)

    tr_auc, W = fit_ridge_lambda(Xq, Y, q_rows[split_tr], y_rows[split_tr],
                                 row_negs, st["query_genes"], cfgA["ridge_lambda_grid"])
    S = (Xq @ W) @ Y.T
    te_uq = np.unique(q_rows[split_te])
    pos = S[q_rows[split_te], y_rows[split_te]]
    neg = np.concatenate([S[qq, row_negs[st["query_genes"][qq]]] for qq in te_uq])
    a, p = auc_p(pos, neg)
    r = {"cell": label, "ridge_train_auc": float(tr_auc), "test_auc": float(a),
         "test_p": float(p), "n_pairs_test": int(split_te.sum()),
         "dim_x": int(Xq.shape[1]), "dim_y": int(Y.shape[1]),
         "n_pool_cell": len(orf_order)}
    print(f"  cell[{label}]: test AUC={a:.4f} p={p:.2g}", flush=True)
    return r


# ------------------------------------------------------------------ track B
def track_b(st, out, W7):
    cfgB = st["cfg"]["trackB"]
    root, cfg = st["root"], st["cfg"]
    CLASSES = cfgB["label_universe"]
    CLS_IDX = {c: i for i, c in enumerate(CLASSES)}
    nC = len(CLASSES)
    sym2row = st["sym2row_h"]

    meta_n = pd.read_csv(root / cfg["inputs"]["higher_metadata"]["path"], sep="\t")
    meta_n = meta_n[meta_n.benchmark_eligible.astype(str).str.lower() == "true"].reset_index(drop=True)
    assert len(meta_n) == 207, len(meta_n)
    genes_n, mat_n = load_matrix(root / cfg["inputs"]["higher_matrix"]["path"])

    # 查询蛋白表示（冻结规则: 成员基因 ESM2 均值 -> l2n; B2 轴再过注入投影）
    Q_es_rows, dropped = [], []
    for _, r in meta_n.iterrows():
        members = [x.strip() for x in re.split(r"[+;,|]", str(r.perturbation_genes)) if x.strip()]
        rows = [sym2row[m] for m in members if m in sym2row]
        if not rows:
            dropped.append(str(r.signature_id))
            Q_es_rows.append(None)
        else:
            Q_es_rows.append(st["X_es_all"][rows].mean(0))
    keep_q = np.array([q is not None for q in Q_es_rows])
    Q_es = l2n(np.stack([q for q in Q_es_rows if q is not None]))
    Q_b2 = l2n(Q_es @ st["W_inj"].T + st["b_inj"])
    meta_n = meta_n[keep_q].reset_index(drop=True)
    n = len(meta_n)
    print(f"  B queries: {n}/207 (dropped no-embedding: {dropped})", flush=True)

    lab_n = meta_n.function_label.values
    grp_n = meta_n.split_group.values
    in_cls = np.array([str(l) in CLS_IDX for l in lab_n])

    # 酵母参考（V5.1 yeast_ref_mask 同配方 ∩ 池覆盖）
    meta_y = st["meta_y"]
    pool_pos = st["pool_pos"]
    ref_idx = [i for i, (_, m) in enumerate(meta_y.iterrows())
               if str(m.split).lower() == "train" and str(m.function_label) in CLS_IDX
               and str(m.target_stable_id) in pool_pos]
    ref_orfs = [str(meta_y.iloc[i].target_stable_id) for i in ref_idx]
    lab_ref = np.array([CLS_IDX[str(meta_y.iloc[i].function_label)] for i in ref_idx])
    ref_rows = [pool_pos[o] for o in ref_orfs]
    Yrb_ref = st["Y_rb"][ref_rows]
    Yes_ref = st["Y_es"][ref_rows]
    print(f"  B yeast refs: {len(ref_idx)}", flush=True)

    # goSlim（lab_gs0 成员, V5.1 逐字机制）
    edges_gs, shared_gs = build_projection_maps(root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)
    Nn = np.stack([mat_n[s] for s in meta_n.signature_id])
    V_n = project_v2_matrix(Nn, genes_n, "human", edges_gs, N_MOD)
    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    Yg_all = np.stack([mat_y[s] for s in meta_y.signature_id])
    V_y_gs = project_v2_matrix(Yg_all, genes_y, "yeast", edges_gs, N_MOD)
    Yg0 = V_y_gs[ref_idx]
    del Nn, Yg_all, mat_y, mat_n

    # 五成员
    S_h = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = Q_b2[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h[i] = cls_scores(Q_b2[i], Href, lab_h, nC)
    S_y_b2 = np.stack([cls_scores(Q_b2[i], Yrb_ref, lab_ref, nC) for i in range(n)])
    S_jrepr = np.stack([cls_scores(Q_b2[i] @ W7, Yrb_ref, lab_ref, nC) for i in range(n)])
    S_y_es = np.stack([cls_scores(Q_es[i], Yes_ref, lab_ref, nC) for i in range(n)])
    lam_lab = cfgB["label_centroid_ridge_lambda"]
    S_lab = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        trm = (grp_n != g) & in_cls
        h_lab = np.array([CLS_IDX[str(x)] for x in lab_n[trm]])
        Ch = label_centroids(V_n[trm], h_lab, nC)
        Cy = label_centroids(Yg0, lab_ref, nC)
        Wl = np.linalg.solve(Ch.T @ Ch + lam_lab * np.eye(Ch.shape[1]), Ch.T @ Cy)
        for i in np.flatnonzero(grp_n == g):
            S_lab[i] = cls_scores(V_n[i] @ Wl, Yg0, lab_ref, nC)

    MEM = {"h_b2": S_h, "y_b2": S_y_b2, "jrepr_b2": S_jrepr,
           "lab_gs0": S_lab, "y_es": S_y_es}
    families = {k: list(v) for k, v in cfgB["families"].items()}
    k_grid = list(cfgB["rrf_k_grid"])

    def rrf(members, k):
        out_ = np.zeros((n, nC))
        for m in members:
            S = MEM[m]
            for i in range(n):
                out_[i] += 1.0 / (k + np.argsort(np.argsort(-S[i])))
        return out_

    rr_cfg = {}
    for fam, members in families.items():
        for k in k_grid:
            Sf = rrf(members, k)
            rr_cfg[(fam, k)] = np.array([
                1.0 / r if (r := rank_of_truth(Sf[i], str(lab_n[i]), CLS_IDX)) > 0 else 0.0
                for i in range(n)])
    rr_single = {m: np.array([
        1.0 / r if (r := rank_of_truth(MEM[m][i], str(lab_n[i]), CLS_IDX)) > 0 else 0.0
        for i in range(n)]) for m in MEM}

    def nested(cand_fams):
        rr_nest = np.full(n, np.nan)
        picks = {}
        for g in pd.unique(grp_n):
            te = grp_n == g
            tr = np.flatnonzero((~te) & in_cls)
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

    rr_fusion, picks_f = nested(list(families))
    rr_comp, picks_c = nested([f for f in families if f.startswith("single:")])
    d = rr_fusion - rr_comp
    p = signflip_p(d[~np.isnan(d)], cfgB["seed_signflip"], cfgB["n_signflip"])
    gate = bool(p < 0.05)
    resB = {"n_queries": n, "dropped_queries": dropped,
            "fusion_mrr": float(np.nanmean(rr_fusion)),
            "comparator_best_single_mrr": float(np.nanmean(rr_comp)),
            "delta": float(np.nanmean(d)), "signflip_p": p, "gate_pass": gate,
            "single_arm_mrr": {m: float(v.mean()) for m, v in rr_single.items()},
            "fusion_picks": picks_f, "comparator_picks": picks_c,
            "v51_reference": {"nested_joint_mrr": 0.348, "signflip_p": 0.00244}}
    print(f"  B fusion MRR={resB['fusion_mrr']:.3f} comp={resB['comparator_best_single_mrr']:.3f} "
          f"Δ={resB['delta']:+.4f} p={p:.4g} gate={gate}", flush=True)
    pd.DataFrame({"signature_id": meta_n.signature_id, "label": lab_n,
                  "split_group": grp_n, "rr_comparator_best_single": rr_comp,
                  "rr_nested_fusion": rr_fusion}).to_csv(
        out / "per_query_rr_trackB.tsv", sep="\t", index=False)
    return resB


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    st = v7_setup(args.config)
    res = {"protocol": st["cfg"]["protocol"],
           "config_sha256": sha256_file(args.config),
           "script_sha256": sha256_file(str(Path(__file__).resolve()))}
    (out / "INPUT_MANIFEST.json").write_text(json.dumps(st["manifest"], indent=1))

    # 主轴 B2（项目门）; 参照轴 esm2_mean（并列报告, 不设项目门）
    resA_b2, W7 = track_a_armset(st, st["X_b2_all"], st["X_b21_all"],
                                 st["Y_rb"], st["Y_rb1"], "b2_primary", out)
    resA_es, _ = track_a_armset(st, st["X_es_all"], st["X_es1_all"],
                                st["Y_es"], st["Y_es1"], "esm2_reference", out)
    res["trackA_primary_b2"] = resA_b2
    res["trackA_reference_esm2"] = resA_es

    # 扩展对集（描述性）: 主轴 train ridge 应用到非 753 的可用 edges
    if st["ext"]:
        hrow = {g: i for i, g in enumerate(st["query_genes"])}
        ext_te = [(h, y) for h, y, s in st["ext"] if s == "test" and h in hrow]
        if len(ext_te) >= 20:
            Xq_b2 = st["X_b2_all"][[st["sym2row_h"][g] for g in st["query_genes"]]]
            S = (Xq_b2 @ W7) @ st["Y_rb"].T
            pos = np.array([S[hrow[h], st["pool_pos"][y]] for h, y in ext_te])
            neg = np.concatenate([S[hrow[h], st["row_negs"][h]] for h, _ in ext_te])
            a, p = auc_p(pos, neg)
            res["trackA_extended_descriptive"] = {"n_test_pairs": len(ext_te),
                                                 "ridge_test_auc": float(a), "p": float(p)}
            print(f"  A[extended, descriptive]: n={len(ext_te)} AUC={a:.4f} p={p:.2g}", flush=True)

    # 附加格（仅对齐器, repr_grid 式）
    res["trackA_extra_cells"] = {
        "b2proj_x_scyeast": aligner_cell(st, st["X_b2_all"], st["Y_scy"], st["scy_pool"],
                                         "X_b2->scyeast_prior"),
        "b2proj_x_yeastESM2": aligner_cell(st, st["X_b2_all"], st["Y_es"], st["pool"],
                                           "X_b2->yeast_esm2"),
        "esm2_x_routeB": aligner_cell(st, st["X_es_all"], st["Y_rb"], st["pool"],
                                      "X_es->route_b"),
        "esm2_x_scyeast": aligner_cell(st, st["X_es_all"], st["Y_scy"], st["scy_pool"],
                                       "X_es->scyeast_prior"),
    }

    res["trackB"] = track_b(st, out, W7)

    (out / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("saved:", out / "result.json")
    print(f"GATE A (B2 primary): primary={resA_b2['gate_primary']} "
          f"secondary={resA_b2['gate_secondary']} | GATE B: {res['trackB']['gate_pass']} "
          f"| ref-axis A (esm2): primary={resA_es['gate_primary']}")


if __name__ == "__main__":
    main()
