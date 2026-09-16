#!/usr/bin/env python3
"""v8_esm2_formal_run — V8 ESM2 轴联合增益正式运行（三臂两门，优化联合臂）。

协议：docs/cross_species_match/V8_ESM2_JOINT_FORMAL_PROTOCOL.md
配置：configs/v8_esm2_joint_formal.json（输入 SHA-256 冻结，运行前强制核对）
机械复用：v7_setup（池/对集/负例/哈希核对）与 V5.1 统计机器。

用法：
  python scripts/v8_esm2_formal_run.py --config configs/v8_esm2_joint_formal.json \
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
from v7_b2_formal_run import (  # noqa: E402
    cls_scores, label_centroids, pct_rows, rank_of_truth, v7_setup)


# ------------------------------------------------------------------ PC 工具
def fit_pc(M):
    mu = M.mean(0)
    _, _, Vt = np.linalg.svd(M - mu, full_matrices=False)
    return mu, Vt


def apply_pc(X, mu, Vt, npc):
    Xc = X - mu
    for i in range(npc):
        Xc = Xc - np.outer(Xc @ Vt[i], Vt[i])
    return l2n(Xc)


# ------------------------------------------------------------------ C 臂配方
def c_recipe(Xq, Y, fx, fy, q_rows, y_rows, row_negs, query_genes, rcp,
             return_W=False):
    """V8 冻结配方: npc 网格 x {ridge, anchored, PLS}, 全部按给定配对的
    train AUC 选参; C = 三成员行百分位秩平均。
    Xq 行序 = query_genes 序 (PC 基 fx/fy 在完整冻结矩阵上拟合后逐行应用,
    与先全矩阵应用再切行等价); Y 行序 = 池序; q_rows/y_rows 为配对索引。
    返回 (S_C, info[, W8, npc*])，S_C 行序 = query_genes 序。"""
    lam_g = rcp["ridge_lambda_grid"]
    anch_g = rcp["anchored_lambda_grid"]
    pls_g = rcp["pls_k_grid"]
    uniq_q = np.unique(q_rows)

    def train_auc_of(S):
        pos = S[q_rows, y_rows]
        neg = np.concatenate([S[qq, row_negs[query_genes[qq]]] for qq in uniq_q])
        return auc_p(pos, neg)[0]

    per_npc = {}
    for npc in rcp["npc_grid"]:
        X1 = apply_pc(Xq, fx[0], fx[1], npc)
        Y1 = apply_pc(Y, fy[0], fy[1], npc)
        Xtr, Ytr = X1[q_rows], Y1[y_rows]
        XtX, XtY = Xtr.T @ Xtr, Xtr.T @ Ytr
        I = np.eye(X1.shape[1])
        cands = {}
        for lam in lam_g:
            W = np.linalg.solve(XtX + lam * I, XtY)
            cands[("ridge", lam)] = (train_auc_of((X1 @ W) @ Y1.T), W, None)
        for lam in anch_g:
            W = np.linalg.solve(XtX + lam * I, XtY + lam * I)
            cands[("anchored", lam)] = (train_auc_of((X1 @ W) @ Y1.T), W, None)
        mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
        my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
        U, _, Vt = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                 full_matrices=False)
        for k in pls_g:
            Px = l2n(((X1 - mx_) / sx_) @ U[:, :k])
            Py = l2n(((Y1 - my_) / sy_) @ Vt[:k].T)
            cands[("pls", k)] = (train_auc_of(Px @ Py.T), None, (mx_, sx_, my_, sy_, U, Vt, k))
        per_npc[npc] = (X1, Y1, cands)

    # npc* = 最优 ridge 的 train AUC 最高的 npc（冻结规则）
    best_npc = max(rcp["npc_grid"],
                   key=lambda n: max(a for (m, _), (a, _, _) in per_npc[n][2].items()
                                     if m == "ridge"))
    X1, Y1, cands = per_npc[best_npc]
    picks = {}
    S_parts = []
    W8 = None
    for mtype in ("ridge", "anchored", "pls"):
        best_key = max((k_ for k_ in cands if k_[0] == mtype), key=lambda k_: cands[k_][0])
        a, W, pls_pack = cands[best_key]
        picks[mtype] = {"hyper": best_key[1], "train_auc": float(a)}
        if W is not None:
            S_parts.append((X1 @ W) @ Y1.T)
            if mtype == "ridge":
                W8 = W
        else:
            mx_, sx_, my_, sy_, U, Vt, k = pls_pack
            Px = l2n(((X1 - mx_) / sx_) @ U[:, :k])
            Py = l2n(((Y1 - my_) / sy_) @ Vt[:k].T)
            S_parts.append(Px @ Py.T)
    S_C = np.mean([pct_rows(S) for S in S_parts], axis=0)
    info = {"npc_selected": int(best_npc), "picks": picks,
            "npc_grid_ridge_train_auc": {
                str(n): float(max(a for (m, _), (a, _, _) in per_npc[n][2].items()
                                  if m == "ridge"))
                for n in rcp["npc_grid"]}}
    if return_W:
        return S_C, info, W8, best_npc
    return S_C, info


# ------------------------------------------------------------------ track A
def track_a(st, out):
    cfgA = st["cfg"]["trackA"]
    rcp = cfgA["armC_recipe"]
    prim = st["prim"]
    q_pos, pool_pos = st["q_pos"], st["pool_pos"]
    row_negs, query_genes = st["row_negs"], st["query_genes"]
    n_pool = st["Y_es"].shape[0]

    Xall = st["X_es_all"]           # l2n ESM2, 810 行 (冻结资产)
    Y = st["Y_es"]                  # l2n 池行
    fx = fit_pc(Xall)
    fy = fit_pc(Y)
    # 一致性自检: apply_pc(npc=1) == remove_pcs(...,1)
    assert np.allclose(apply_pc(Xall, *fx, 1), remove_pcs(Xall, 1), atol=1e-10)

    hrow = [st["sym2row_h"][g] for g in query_genes]
    Xq = Xall[hrow]
    q_rows = np.array([q_pos[h] for h, _, _, _ in prim])
    y_rows = np.array([pool_pos[y] for _, y, _, _ in prim])
    split_tr = np.array([s == "train" for _, _, s, _ in prim])
    split_te = np.array([s == "test" for _, _, s, _ in prim])
    folds = np.array([f for _, _, _, f in prim])

    S_A = Xq @ Y.T
    S_B = apply_pc(Xq, *fx, 1) @ apply_pc(Y, *fy, 1).T
    S_C, info, W8, npc_star = c_recipe(
        Xq, Y, fx, fy, q_rows[split_tr], y_rows[split_tr], row_negs,
        query_genes, rcp, return_W=True)
    resA = {"repr": "esm2_axis_v8", "n_pairs_primary": len(prim), "n_pool": n_pool,
            "armC_recipe_info": info}
    print(f"  A recipe: npc*={info['npc_selected']} picks={info['picks']}", flush=True)
    te_uq = np.unique(q_rows[split_te])

    def arm_eval(S, name):
        pos = S[q_rows[split_te], y_rows[split_te]]
        neg = np.concatenate([S[qq, row_negs[query_genes[qq]]] for qq in te_uq])
        a, p = auc_p(pos, neg)
        resA[f"{name}_test_auc"] = a
        resA[f"{name}_test_p"] = p
        print(f"  A {name}: test auc={a:.4f} p={p:.2g}", flush=True)

    # 三臂分数矩阵统一为 query_genes 行序
    S_Cq = S_C
    for nm, S in (("armA_no_align", S_A), ("armB_unsup_pc1", S_B), ("armC_joint", S_Cq)):
        arm_eval(S, nm)

    def boot_delta(Sa, Sb, seed):
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

    for nm, Sa, Sb in (("C_vs_A", S_Cq, S_A), ("C_vs_B", S_Cq, S_B), ("B_vs_A", S_B, S_A)):
        m, lo, hi = boot_delta(Sa, Sb, cfgA["seed_bootstrap"])
        resA[f"boot_{nm}"] = {"mean": m, "ci95_low": lo, "ci95_high": hi}
        print(f"  A ΔAUC {nm}: {m:+.4f} CI[{lo:+.4f},{hi:+.4f}]", flush=True)

    # OOF（冻结 fold; 每折内完整重跑配方选择）
    oof = {nm: np.full(len(prim), np.nan) for nm in ("A", "B", "C")}
    oof_hit = {nm: np.full(len(prim), np.nan) for nm in ("B", "C")}

    def pct_ranks(S, mask):
        return np.array([(np.sum(S[qq] > S[qq, yy]) + 1) / n_pool
                         for qq, yy in zip(q_rows[mask], y_rows[mask])])

    for f in sorted(set(folds.tolist())):
        m_tr, m_te = folds != f, folds == f
        S_Cf, _ = c_recipe(Xq, Y, fx, fy, q_rows[m_tr], y_rows[m_tr],
                           row_negs, query_genes, rcp)
        for nm, S in (("A", S_A), ("B", S_B), ("C", S_Cf)):
            oof[nm][m_te] = pct_ranks(S, m_te)
        for i in np.flatnonzero(m_te):
            g = query_genes[q_rows[i]]
            partners = [pool_pos[y2] for y2 in st["h_partner"].get(g, ()) if y2 in pool_pos]
            if not partners:
                continue
            oof_hit["B"][i] = min(int(np.sum(S_B[q_rows[i]] > S_B[q_rows[i], c_])) + 1
                                  for c_ in partners) <= 10
            oof_hit["C"][i] = min(int(np.sum(S_Cf[q_rows[i]] > S_Cf[q_rows[i], c_])) + 1
                                  for c_ in partners) <= 10
    valid = ~np.isnan(oof["A"])
    sf_seed = cfgA["seed_signflip_base"]
    for nm_a, nm_b in (("C", "A"), ("C", "B")):
        d = oof[nm_b][valid] - oof[nm_a][valid]
        sf_seed += 1
        p_ = signflip_p(d, sf_seed, cfgA["n_signflip"])
        resA[f"oof_{nm_a}_vs_{nm_b}"] = {"gain": float(d.mean()), "signflip_p": p_}
        print(f"  A OOF {nm_a}-vs-{nm_b}: {d.mean():+.4f} p={p_:.3g}", flush=True)
    hb = oof_hit["B"][~np.isnan(oof_hit["B"])]
    hc = oof_hit["C"][~np.isnan(oof_hit["C"])]
    d_hit = hc.astype(float) - hb.astype(float)
    sf_seed += 1
    resA["any_hit10"] = {"rate_B": float(hb.mean()), "rate_C": float(hc.mean()),
                         "signflip_p": signflip_p(d_hit, sf_seed, cfgA["n_signflip"])}
    print(f"  A any-hit@10: B={hb.mean():.4f} C={hc.mean():.4f} "
          f"p={resA['any_hit10']['signflip_p']:.3g}", flush=True)

    resA["gate_primary"] = bool(resA["armC_joint_test_p"] < 0.05 and
                                resA["boot_C_vs_A"]["ci95_low"] > 0)
    resA["gate_secondary"] = bool(resA["oof_C_vs_A"]["signflip_p"] < 0.05 or
                                  resA["any_hit10"]["signflip_p"] < 0.05)
    up = (resA["boot_C_vs_B"]["ci95_low"] > 0 and
          resA["oof_C_vs_B"]["signflip_p"] < 0.05)
    resA["upgraded_claim_C_beats_B"] = bool(up)
    print(f"  A gates: primary={resA['gate_primary']} secondary={resA['gate_secondary']} "
          f"| upgraded C>B claim: {up}", flush=True)

    pd.DataFrame({"human": [h for h, _, _, _ in prim],
                  "yeast_orf": [y for _, y, _, _ in prim],
                  "split": [s for _, _, s, _ in prim], "fold": folds,
                  "oof_pct_rank_A": oof["A"], "oof_pct_rank_B": oof["B"],
                  "oof_pct_rank_C": oof["C"]}).to_csv(
        out / "per_query_ranks_trackA.tsv", sep="\t", index=False)
    return resA, {"W8": W8, "npc_star": npc_star, "fx": fx, "fy": fy}


# ------------------------------------------------------------------ track B
def track_b(st, out, ta_state):
    cfgB = st["cfg"]["trackB"]
    root, cfg = st["root"], st["cfg"]
    CLASSES = cfgB["label_universe"]
    CLS_IDX = {c: i for i, c in enumerate(CLASSES)}
    nC = len(CLASSES)
    sym2row = st["sym2row_h"]
    W8, npc_star, fx, fy = ta_state["W8"], ta_state["npc_star"], ta_state["fx"], ta_state["fy"]

    meta_n = pd.read_csv(root / cfg["inputs"]["higher_metadata"]["path"], sep="\t")
    meta_n = meta_n[meta_n.benchmark_eligible.astype(str).str.lower() == "true"].reset_index(drop=True)
    assert len(meta_n) == 207, len(meta_n)
    genes_n, mat_n = load_matrix(root / cfg["inputs"]["higher_matrix"]["path"])

    Q_raw_rows, dropped = [], []
    for _, r in meta_n.iterrows():
        members = [x.strip() for x in re.split(r"[+;,|]", str(r.perturbation_genes)) if x.strip()]
        rows = [sym2row[m] for m in members if m in sym2row]
        if not rows:
            dropped.append(str(r.signature_id))
            Q_raw_rows.append(None)
        else:
            Q_raw_rows.append(st["X_es_all"][rows].mean(0))
    keep_q = np.array([q is not None for q in Q_raw_rows])
    Q_raw = l2n(np.stack([q for q in Q_raw_rows if q is not None]))
    Q1 = apply_pc(Q_raw, fx[0], fx[1], 1)
    QN = apply_pc(Q_raw, fx[0], fx[1], npc_star)
    Q_b2 = l2n(Q_raw @ st["W_inj"].T + st["b_inj"])
    meta_n = meta_n[keep_q].reset_index(drop=True)
    n = len(meta_n)
    print(f"  B queries: {n}/207 (dropped: {dropped})", flush=True)

    lab_n = meta_n.function_label.values
    grp_n = meta_n.split_group.values
    in_cls = np.array([str(l) in CLS_IDX for l in lab_n])

    meta_y = st["meta_y"]
    pool_pos = st["pool_pos"]
    ref_idx = [i for i, (_, m) in enumerate(meta_y.iterrows())
               if str(m.split).lower() == "train" and str(m.function_label) in CLS_IDX
               and str(m.target_stable_id) in pool_pos]
    ref_rows = [pool_pos[str(meta_y.iloc[i].target_stable_id)] for i in ref_idx]
    lab_ref = np.array([CLS_IDX[str(meta_y.iloc[i].function_label)] for i in ref_idx])
    Y1_ref = apply_pc(st["Y_es"], *fy, 1)[ref_rows]
    YN_ref = apply_pc(st["Y_es"], *fy, npc_star)[ref_rows]
    Yrb_ref = st["Y_rb"][ref_rows]

    # lab_gs0（V5.1 逐字机制）
    edges_gs, shared_gs = build_projection_maps(root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)
    Nn = np.stack([mat_n[s] for s in meta_n.signature_id])
    V_n = project_v2_matrix(Nn, genes_n, "human", edges_gs, N_MOD)
    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    Yg_all = np.stack([mat_y[s] for s in meta_y.signature_id])
    V_y_gs = project_v2_matrix(Yg_all, genes_y, "yeast", edges_gs, N_MOD)
    Yg0 = V_y_gs[ref_idx]
    del Nn, Yg_all, mat_y, mat_n

    S_h = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = Q1[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h[i] = cls_scores(Q1[i], Href, lab_h, nC)
    S_y1 = np.stack([cls_scores(Q1[i], Y1_ref, lab_ref, nC) for i in range(n)])
    S_j8 = np.stack([cls_scores(QN[i] @ W8, YN_ref, lab_ref, nC) for i in range(n)])
    S_b2 = np.stack([cls_scores(Q_b2[i], Yrb_ref, lab_ref, nC) for i in range(n)])
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

    MEM = {"h_es1": S_h, "y_es1": S_y1, "jrepr8": S_j8, "lab_gs0": S_lab, "y_b2": S_b2}
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
            "v7_reference": {"fusion": 0.369, "comparator": 0.351, "p": 0.124}}
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

    resA, ta_state = track_a(st, out)
    res["trackA"] = resA
    res["trackB"] = track_b(st, out, ta_state)

    (out / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("saved:", out / "result.json")
    print(f"GATE A: primary={resA['gate_primary']} secondary={resA['gate_secondary']} "
          f"| upgraded C>B: {resA['upgraded_claim_C_beats_B']} "
          f"| GATE B: {res['trackB']['gate_pass']}")


if __name__ == "__main__":
    main()
