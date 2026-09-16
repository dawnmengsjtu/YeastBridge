#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EV1a — 基本结构外样本（Replogle held-out 基因查询，Norman 标签直传）。
V9 三臂五成员机器 verbatim；LOGO 组 = 单基因；判据见 EV_PROTOCOL.md §EV1。"""
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

YB = Path("/public/home/mengxl/dzy/YeastBridge_0912")
EV = Path("/public/home/mengxl/dzy/YeastBridge_0912")
OUT = EV / "results_a6000/externalvalidation/ev1a"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(YB / "scripts"))
from v5_joint_formal_run import (build_projection_maps, l2n,  # noqa: E402
                                 project_v2_matrix, sha256_file, signflip_p)
from v7_b2_formal_run import cls_scores, label_centroids, pct_rows, v7_setup  # noqa: E402
from v8_esm2_formal_run import apply_pc, fit_pc  # noqa: E402
from cross_species_match import load_matrix  # noqa: E402

NPC, LAM_R, LAM_A, PLS = 3, 1.0, 1.0, 8


def rank_of_truth(scores, truth, CI):
    if truth not in CI:
        return 0
    return int(np.where(np.argsort(-scores, kind="stable") == CI[truth])[0][0]) + 1


def main():
    v9cfg = json.loads((YB / "configs_a6000_frozen/v9_esm2_basic.json").read_text())
    cfgB = v9cfg["trackB"]
    CLASSES = cfgB["label_universe"]
    CI = {c: i for i, c in enumerate(CLASSES)}
    nC = len(CLASSES)
    # 冻结超参断言（与 EJ 方法冻结一致）
    ej = json.loads((YB / "configs_a6000_frozen/esm2_joint_tasks.json")
                    .read_text())["method"]
    assert (ej["npc_aligned"], ej["ridge_lambda"], ej["anchored_lambda"],
            ej["pls_k"]) == (NPC, LAM_R, LAM_A, PLS)

    st = v7_setup(str(YB / "configs/ev_v7_mirror.json"))
    sym2row, pool_pos = st["sym2row_h"], st["pool_pos"]
    Xall, Y = st["X_es_all"], st["Y_es"]
    fx, fy = fit_pc(Xall), fit_pc(Y)
    mir = YB / "raw/externalvalidation/frozen_inputs"

    # ---- 查询嵌入（810 ∪ EV0x 补推理表）----
    uni = pd.read_csv(EV / "raw/externalvalidation/ev0x_embed_merged/index.tsv", sep="\t").fillna("")
    U = np.load(EV / "raw/externalvalidation/ev0x_embed_merged/esm2_mean_fp32.npy").astype(np.float64)
    urow = {str(s).upper(): i for i, s in enumerate(uni["common"])}
    q = pd.read_csv(OUT / "ev1a_query_set.tsv", sep="\t")
    rows = []
    for g in q.gene:
        gu = str(g).upper()
        if gu in sym2row:
            rows.append(Xall[sym2row[gu]])
        else:
            rows.append(l2n(U[[urow[gu]]])[0])
    Q_raw = l2n(np.stack(rows))
    n = len(q)
    lab = q.function_label.astype(str).values
    grp = q.gene.astype(str).values
    in_cls = np.array([l in CI for l in lab])
    Q1 = apply_pc(Q_raw, fx[0], fx[1], 1)
    QN = apply_pc(Q_raw, fx[0], fx[1], NPC)
    Q_b2 = l2n(Q_raw @ st["W_inj"].T + st["b_inj"])

    # ---- 酵母带标签参考 ----
    meta_y, prim = st["meta_y"], st["prim"]
    ref_idx = [i for i, (_, m) in enumerate(meta_y.iterrows())
               if str(m.split).lower() == "train" and str(m.function_label) in CI
               and str(m.target_stable_id) in pool_pos]
    ref_rows = [pool_pos[str(meta_y.iloc[i].target_stable_id)] for i in ref_idx]
    lab_ref = np.array([CI[str(meta_y.iloc[i].function_label)] for i in ref_idx])
    YN_full = apply_pc(Y, fy[0], fy[1], NPC)
    Y1 = apply_pc(Y, fy[0], fy[1], 1)[ref_rows]
    YN = YN_full[ref_rows]
    Yrb = st["Y_rb"][ref_rows]

    # ---- j9（冻结超参，官方 train 对确定性重拟合）----
    q_pos = st["q_pos"]
    hrow = [st["sym2row_h"][g] for g in st["query_genes"]]
    X753 = apply_pc(Xall[hrow], fx[0], fx[1], NPC)
    qr = np.array([q_pos[h] for h, _, _, _ in prim])
    yr = np.array([pool_pos[y] for _, y, _, _ in prim])
    trm = np.array([s == "train" for _, _, s, _ in prim])
    Xtr, Ytr = X753[qr[trm]], YN_full[yr[trm]]
    XtX, XtY = Xtr.T @ Xtr, Xtr.T @ Ytr
    I = np.eye(Xtr.shape[1])
    W_r = np.linalg.solve(XtX + LAM_R * I, XtY)
    W_a = np.linalg.solve(XtX + LAM_A * I, XtY + LAM_A * I)
    mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    Usv, _, Vt = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                               full_matrices=False)
    S_r = (QN @ W_r) @ YN.T
    S_a = (QN @ W_a) @ YN.T
    S_j9 = np.mean([pct_rows(S_r), pct_rows(S_a),
                    pct_rows(l2n(((QN - mx_) / sx_) @ Usv[:, :PLS]) @
                             l2n(((YN - my_) / sy_) @ Vt[:PLS].T).T)], axis=0)
    S9 = np.full((n, nC), -np.inf)
    for c in range(nC):
        m = lab_ref == c
        if m.any():
            S9[:, c] = S_j9[:, m].max(axis=1)

    # ---- lab_gs0：K562 矩阵查询行投影（先镜像 h5ad，自包含）----
    v5cfg = json.loads((YB / "configs_a6000_frozen/v5_1_joint_formal.json").read_text())
    v5stage = Path("/public/home/mengxl/dzy/YeastBridge_0912/raw/externalvalidation/frozen_inputs")
    rel = v5cfg["inputs"]["k562_h5ad"]["path"]
    h5src = v5stage / rel
    h5dst = mir / rel
    if not h5dst.exists():
        h5dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(h5src, h5dst)
    assert sha256_file(str(h5dst)) == v5cfg["inputs"]["k562_h5ad"]["sha256"]
    import anndata as ad
    a = ad.read_h5ad(h5dst)
    genes_k = list(a.var.gene_name)
    obsn = list(a.obs.index)
    pat = re.compile(r"^\d+_([^_]+)_(P1|P1P2)_")
    keep, tgt = [], []
    for name in obsn:
        m = pat.match(name)
        if m and m.group(2) in {"P1", "P2", "P1P2"}:
            keep.append(name)
            tgt.append(m.group(1))
    Xp = a.X[[obsn.index(x) for x in keep]]
    Xp = Xp.toarray() if hasattr(Xp, "toarray") else np.asarray(Xp)
    Xp = np.nan_to_num(Xp, nan=0.0)
    col_ok = np.abs(Xp).max(0) < 1e30
    Xp = Xp[:, col_ok]
    genes_k = [g for g, ok in zip(genes_k, col_ok) if ok]
    t2row = {}
    for j, t in enumerate(tgt):
        t2row.setdefault(t.upper(), j)
    ridx = [t2row[g.upper()] for g in q.gene if g.upper() in t2row]
    have = [g for g in q.gene.astype(str) if g.upper() in t2row]
    Xq_k = Xp[[t2row[g.upper()] for g in have]]
    gproj = mir / "match_diag_20260908/k562_proj/gene_goslim_projection.tsv.gz"
    edges_gs, shared_gs = build_projection_maps(gproj)
    V_q = project_v2_matrix(Xq_k, genes_k, "human", edges_gs, len(shared_gs))
    genes_y, mat_y = load_matrix(mir / "gse_go_three_arm_v3_open_set/"
                                 "inputs_v3_open_set/yeast_matrix.tsv.gz")
    Yg_all = np.stack([mat_y[s] for s in meta_y.signature_id])
    V_y = project_v2_matrix(Yg_all, genes_y, "yeast", edges_gs, len(shared_gs))
    Yg0 = V_y[ref_idx]
    lam_lab = cfgB["label_centroid_ridge_lambda"]
    # 行对齐：V_q 行序 = have 序；成员分数补缺行用 -inf（与 V9 丢查询不同，此处保形）
    S_lab_all = np.full((n, nC), -np.inf)
    pos_of = {g: i for i, g in enumerate(q.gene.astype(str))}
    for g, vq in zip(have, V_q):
        i = pos_of[g]
        Ch = label_centroids(V_q[[pos_of[x] for x in have if x != g]],
                             np.array([CI[lab[pos_of[x]]] for x in have
                                       if x != g and in_cls[pos_of[x]]]), nC) \
            if False else None
        trm2 = np.array([x != g and in_cls[pos_of[x]] for x in have])
        Ch = label_centroids(V_q[trm2], np.array([CI[lab[pos_of[x]]]
                                                  for x in np.array(have)[trm2]]), nC)
        Cy = label_centroids(Yg0, lab_ref, nC)
        Wl = np.linalg.solve(Ch.T @ Ch + lam_lab * np.eye(Ch.shape[1]), Ch.T @ Cy)
        S_lab_all[i] = cls_scores(vq @ Wl, Yg0, lab_ref, nC)

    # ---- 五成员 + RRF 嵌套 ----
    S_h = np.full((n, nC), -np.inf)
    for g in pd.unique(grp):
        hm = (grp != g) & in_cls
        Href, lh = Q1[hm], np.array([CI[x] for x in lab[hm]])
        for i in np.flatnonzero(grp == g):
            S_h[i] = cls_scores(Q1[i], Href, lh, nC)
    S_y1 = np.stack([cls_scores(Q1[i], Y1, lab_ref, nC) for i in range(n)])
    S_b2 = np.stack([cls_scores(Q_b2[i], Yrb, lab_ref, nC) for i in range(n)])
    MEM = {"h_es1": S_h, "y_es1": S_y1, "j9": S9, "lab_gs0": S_lab_all, "y_b2": S_b2}
    fams = {k: list(v) for k, v in cfgB["families"].items()}
    kg = list(cfgB["rrf_k_grid"])

    def rrf(members, k):
        o = np.zeros((n, nC))
        for m_ in members:
            S = MEM[m_]
            for i in range(n):
                o[i] += 1.0 / (k + np.argsort(np.argsort(-S[i])))
        return o

    rr_cfg = {}
    for f, ms in fams.items():
        for k in kg:
            Sf = rrf(ms, k)
            rr_cfg[(f, k)] = np.array([
                1.0 / r if (r := rank_of_truth(Sf[i], lab[i], CI)) > 0 else 0.0
                for i in range(n)])

    def nested(cands):
        rr = np.full(n, np.nan)
        picks = {}
        for g in pd.unique(grp):
            te = grp == g
            trn = np.flatnonzero((~te) & in_cls)
            sc = []
            for f in cands:
                for k in kg:
                    sc.append(((-float(rr_cfg[(f, k)][trn].mean()), len(fams[f]), k),
                               (f, k)))
            sc.sort(key=lambda x: x[0])
            f_, k_ = sc[0][1]
            picks[f"{f_}|k={k_}"] = picks.get(f"{f_}|k={k_}", 0) + 1
            for i in np.flatnonzero(te):
                rr[i] = rr_cfg[(f_, k_)][i]
        return rr, picks

    rr_j, pj = nested(list(fams))
    rr_c, pc_ = nested([f for f in fams if f.startswith("single:")])
    d1 = rr_j - rr_c
    p1 = signflip_p(d1[~np.isnan(d1)], cfgB["seed_signflip_base"] + 1,
                    cfgB["n_signflip"])
    single = {m: np.array([1.0 / r if (r := rank_of_truth(MEM[m][i], lab[i], CI)) > 0
                           else 0.0 for i in range(n)]) for m in MEM}
    dh, dy = rr_j - single["h_es1"], rr_j - single["y_es1"]
    ph = signflip_p(dh[~np.isnan(dh)], cfgB["seed_signflip_base"] + 2,
                    cfgB["n_signflip"])
    py_ = signflip_p(dy[~np.isnan(dy)], cfgB["seed_signflip_base"] + 3,
                     cfgB["n_signflip"])
    jm, hm_, ym = (float(np.nanmean(rr_j)), float(single["h_es1"].mean()),
                   float(single["y_es1"].mean()))
    res = {"mode": "ev1a_basic_structure", "n_queries": n,
           "human_only_mrr": hm_, "yeast_only_mrr": ym,
           "joint_nested_mrr": jm,
           "comparator_best_single_mrr": float(np.nanmean(rr_c)),
           "single_arm_mrr": {m: float(v.mean()) for m, v in single.items()},
           "G1_joint_vs_best_single": {"delta": float(np.nanmean(d1)),
                                       "signflip_p": p1,
                                       "pass": bool(p1 < 0.05)},
           "G2_joint_significantly_best": {
               "delta_vs_human": float(dh.mean()), "p_vs_human": ph,
               "delta_vs_yeast": float(dy.mean()), "p_vs_yeast": py_,
               "pass": bool(ph < 0.05 and py_ < 0.05 and dh.mean() > 0
                            and dy.mean() > 0)},
           "joint_picks": pj,
           "k562_rows_found": len(have)}
    pd.DataFrame({"gene": q.gene, "label": lab,
                  "rr_human_only": single["h_es1"], "rr_yeast_only": single["y_es1"],
                  "rr_comparator_best_single": rr_c,
                  "rr_joint_nested": rr_j}).to_csv(
        OUT / "per_query_rr_ev1a.tsv", sep="	", index=False)
    (OUT / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps({k: v for k, v in res.items()
                      if k not in ("single_arm_mrr", "joint_picks")},
                     ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
