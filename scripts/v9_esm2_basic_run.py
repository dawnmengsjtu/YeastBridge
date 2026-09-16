#!/usr/bin/env python3
"""v9_esm2_basic_run — V9 esm2 轴基本结构三臂正式运行（human-only / yeast-only / joint）。

用户澄清（2026-09-12 深夜）：三臂两门 = v3 基本结构语义（README 的
human_only / yeast_only / joint），交付标准 = 一个方法，joint 显著最优，
且跨物种匹配显著。V9 把两件事合并为单一记录：

  G3 跨物种匹配显著：检索三臂（V8 冻结配方逐字重跑，含门 A 主/次门）；
  G1/G2 joint 显著最优：Norman 19 类基本结构三臂——
    human-only = h_es1（LOGO 人侧参考）
    yeast-only = y_es1（酵母带标签参考）
    joint      = 嵌套 (family,k) 选择的 RRF 融合，家族集合首次包含中间家族
                 pair_hy / triple_hyj9 / quad_hyjb（V8 只测过单臂与 all5，
                 all5 被弱成员稀释——动机来自 V8 已登记单臂 MRR，事后披露）。
    j9 = V8 过门 A 的完整 C 配方算子（npc* 空间 ridge + identity-anchored
         ridge + PLS 秩平均，strain 级分数后按类取 max），替代 V8 的单 ridge
         jrepr8，作为 joint 中的配对监督跨物种成员。

协议：docs/cross_species_match/V9_ESM2_BASIC_STRUCTURE_PROTOCOL.md
配置：configs/v9_esm2_basic.json

用法：
  python scripts/v9_esm2_basic_run.py --config configs/v9_esm2_basic.json \
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
    build_projection_maps, l2n, project_v2_matrix, sha256_file, signflip_p)
from v7_b2_formal_run import (  # noqa: E402
    cls_scores, label_centroids, pct_rows, rank_of_truth, v7_setup)
from v8_esm2_formal_run import apply_pc, fit_pc, track_a  # noqa: E402


def basic_structure(st, out, ta_state, resA):
    cfgB = st["cfg"]["trackB"]
    root, cfg = st["root"], st["cfg"]
    CLASSES = cfgB["label_universe"]
    CLS_IDX = {c: i for i, c in enumerate(CLASSES)}
    nC = len(CLASSES)
    sym2row = st["sym2row_h"]
    npc_star = ta_state["npc_star"]
    fx, fy = ta_state["fx"], ta_state["fy"]
    picks = resA["armC_recipe_info"]["picks"]

    # ---- Norman 查询（V8 逐字口径） ----
    meta_n = pd.read_csv(root / cfg["inputs"]["higher_metadata"]["path"], sep="\t")
    meta_n = meta_n[meta_n.benchmark_eligible.astype(str).str.lower() == "true"].reset_index(drop=True)
    assert len(meta_n) == 207, len(meta_n)
    genes_n, mat_n = load_matrix(root / cfg["inputs"]["higher_matrix"]["path"])
    Q_rows, dropped = [], []
    for _, r in meta_n.iterrows():
        members = [x.strip() for x in re.split(r"[+;,|]", str(r.perturbation_genes)) if x.strip()]
        rows = [sym2row[m] for m in members if m in sym2row]
        if not rows:
            dropped.append(str(r.signature_id))
            Q_rows.append(None)
        else:
            Q_rows.append(st["X_es_all"][rows].mean(0))
    keep_q = np.array([q is not None for q in Q_rows])
    Q_raw = l2n(np.stack([q for q in Q_rows if q is not None]))
    Q1 = apply_pc(Q_raw, fx[0], fx[1], 1)
    QN = apply_pc(Q_raw, fx[0], fx[1], npc_star)
    Q_b2 = l2n(Q_raw @ st["W_inj"].T + st["b_inj"])
    meta_n = meta_n[keep_q].reset_index(drop=True)
    n = len(meta_n)
    lab_n = meta_n.function_label.values
    grp_n = meta_n.split_group.values
    in_cls = np.array([str(l) in CLS_IDX for l in lab_n])
    print(f"  basic-structure queries: {n}/207 (dropped: {dropped})", flush=True)

    # ---- 酵母带标签参考（V5.1 yeast_ref_mask 同配方 ∩ 池） ----
    meta_y = st["meta_y"]
    pool_pos = st["pool_pos"]
    ref_idx = [i for i, (_, m) in enumerate(meta_y.iterrows())
               if str(m.split).lower() == "train" and str(m.function_label) in CLS_IDX
               and str(m.target_stable_id) in pool_pos]
    ref_rows = [pool_pos[str(meta_y.iloc[i].target_stable_id)] for i in ref_idx]
    lab_ref = np.array([CLS_IDX[str(meta_y.iloc[i].function_label)] for i in ref_idx])
    Y1_ref = apply_pc(st["Y_es"], *fy, 1)[ref_rows]
    YN_ref = apply_pc(st["Y_es"], fy[0], fy[1], npc_star)[ref_rows]
    Yrb_ref = st["Y_rb"][ref_rows]
    print(f"  yeast labeled refs: {len(ref_idx)}", flush=True)

    # ---- j9：V8 过门 A 的完整 C 配方算子（同 train 对、同选定超参，确定性重拟合） ----
    prim = st["prim"]
    q_pos = st["q_pos"]
    hrow = [st["sym2row_h"][g] for g in st["query_genes"]]
    Xq1 = apply_pc(st["X_es_all"][hrow], fx[0], fx[1], npc_star)
    Y1m = apply_pc(st["Y_es"], fy[0], fy[1], npc_star)
    q_rows = np.array([q_pos[h] for h, _, _, _ in prim])
    y_rows = np.array([pool_pos[y] for _, y, _, _ in prim])
    split_tr = np.array([s == "train" for _, _, s, _ in prim])
    Xtr, Ytr = Xq1[q_rows[split_tr]], Y1m[y_rows[split_tr]]
    XtX, XtY = Xtr.T @ Xtr, Xtr.T @ Ytr
    I = np.eye(Xq1.shape[1])
    lam_r = float(picks["ridge"]["hyper"])
    lam_a = float(picks["anchored"]["hyper"])
    k_pls = int(picks["pls"]["hyper"])
    W_r = np.linalg.solve(XtX + lam_r * I, XtY)
    W_a = np.linalg.solve(XtX + lam_a * I, XtY + lam_a * I)
    mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    U, _, Vt = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                             full_matrices=False)
    S_r = (QN @ W_r) @ YN_ref.T
    S_a = (QN @ W_a) @ YN_ref.T
    Px = l2n(((QN - mx_) / sx_) @ U[:, :k_pls])
    Py = l2n(((YN_ref - my_) / sy_) @ Vt[:k_pls].T)
    S_p = Px @ Py.T
    S_j9_strain = np.mean([pct_rows(S_r), pct_rows(S_a), pct_rows(S_p)], axis=0)
    S_j9 = np.full((n, nC), -np.inf)
    for c in range(nC):
        m = lab_ref == c
        if m.any():
            S_j9[:, c] = S_j9_strain[:, m].max(axis=1)

    # ---- 其余四成员（V8 逐字机制） ----
    S_h = np.full((n, nC), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = Q1[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h[i] = cls_scores(Q1[i], Href, lab_h, nC)
    S_y1 = np.stack([cls_scores(Q1[i], Y1_ref, lab_ref, nC) for i in range(n)])
    S_b2 = np.stack([cls_scores(Q_b2[i], Yrb_ref, lab_ref, nC) for i in range(n)])
    edges_gs, shared_gs = build_projection_maps(root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)
    Nn = np.stack([mat_n[s] for s in meta_n.signature_id])
    V_n = project_v2_matrix(Nn, genes_n, "human", edges_gs, N_MOD)
    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    Yg_all = np.stack([mat_y[s] for s in meta_y.signature_id])
    V_y_gs = project_v2_matrix(Yg_all, genes_y, "yeast", edges_gs, N_MOD)
    Yg0 = V_y_gs[ref_idx]
    del Nn, Yg_all, mat_y, mat_n
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

    MEM = {"h_es1": S_h, "y_es1": S_y1, "j9": S_j9, "lab_gs0": S_lab, "y_b2": S_b2}
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
        picks_ = {}
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
            picks_[f"{fam_pick}|k={k_pick}"] = picks_.get(f"{fam_pick}|k={k_pick}", 0) + 1
            for i in np.flatnonzero(te):
                rr_nest[i] = rr_cfg[(fam_pick, k_pick)][i]
        return rr_nest, picks_

    rr_joint, picks_j = nested(list(families))
    rr_comp, picks_c = nested([f for f in families if f.startswith("single:")])
    sf = cfgB["seed_signflip_base"]
    n_sf = cfgB["n_signflip"]

    d_comp = rr_joint - rr_comp
    p_g1 = signflip_p(d_comp[~np.isnan(d_comp)], sf + 1, n_sf)
    d_h = rr_joint - rr_single["h_es1"]
    d_y = rr_joint - rr_single["y_es1"]
    p_h = signflip_p(d_h, sf + 2, n_sf)
    p_y = signflip_p(d_y, sf + 3, n_sf)
    jm, hm, ym = (float(np.nanmean(rr_joint)), float(rr_single["h_es1"].mean()),
                  float(rr_single["y_es1"].mean()))
    g1 = bool(p_g1 < 0.05)
    g2 = bool(p_h < 0.05 and p_y < 0.05 and d_h.mean() > 0 and d_y.mean() > 0
              and jm > hm and jm > ym)
    g3 = bool(resA["gate_primary"] and resA["gate_secondary"]
              and resA["armC_joint_test_p"] < 0.05
              and resA["armB_unsup_pc1_test_p"] < 0.05)

    resB = {
        "task": "v3 基本结构三臂 (human-only / yeast-only / joint), esm2 轴",
        "n_queries": n, "dropped_queries": dropped,
        "human_only_mrr": hm, "yeast_only_mrr": ym,
        "joint_nested_mrr": jm,
        "comparator_best_single_mrr": float(np.nanmean(rr_comp)),
        "single_arm_mrr": {m: float(v.mean()) for m, v in rr_single.items()},
        "G1_joint_vs_best_single": {"delta": float(np.nanmean(d_comp)),
                                    "signflip_p": p_g1, "pass": g1},
        "G2_joint_significantly_best": {
            "delta_vs_human_only": float(d_h.mean()), "p_vs_human_only": p_h,
            "delta_vs_yeast_only": float(d_y.mean()), "p_vs_yeast_only": p_y,
            "pass": g2},
        "joint_picks": picks_j, "comparator_picks": picks_c,
        "v8_reference": {"fusion_all5": 0.352, "h_es1": 0.365, "y_es1": 0.356},
    }
    print(f"  basic structure: human-only={hm:.3f} yeast-only={ym:.3f} "
          f"joint={jm:.3f} comp={resB['comparator_best_single_mrr']:.3f}", flush=True)
    print(f"  G1 Δ={resB['G1_joint_vs_best_single']['delta']:+.4f} p={p_g1:.4g} pass={g1}", flush=True)
    print(f"  G2 vs h: Δ={d_h.mean():+.4f} p={p_h:.4g} | vs y: Δ={d_y.mean():+.4f} "
          f"p={p_y:.4g} pass={g2}", flush=True)
    print(f"  joint picks: {picks_j}", flush=True)

    pd.DataFrame({"signature_id": meta_n.signature_id, "label": lab_n,
                  "split_group": grp_n, "rr_human_only": rr_single["h_es1"],
                  "rr_yeast_only": rr_single["y_es1"], "rr_j9": rr_single["j9"],
                  "rr_comparator_best_single": rr_comp,
                  "rr_joint_nested": rr_joint}).to_csv(
        out / "per_query_rr_basic_structure.tsv", sep="\t", index=False)
    return resB, {"G1": g1, "G2": g2, "G3": g3}


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

    # G3: 检索三臂（V8 冻结配方逐字重跑, 兼作复现自检）
    resA, ta_state = track_a(st, out)
    res["trackA_retrieval"] = resA

    # G1/G2: 基本结构三臂
    resB, gates = basic_structure(st, out, ta_state, resA)
    res["trackB_basic_structure"] = resB
    res["gates"] = {"G1_joint_vs_best_single": gates["G1"],
                    "G2_joint_significantly_best": gates["G2"],
                    "G3_cross_species_matching": gates["G3"]}

    (out / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("saved:", out / "result.json")
    print(f"GATES: G1(joint>最优单臂)={gates['G1']} "
          f"G2(joint显著最优)={gates['G2']} G3(跨物种匹配显著)={gates['G3']}")


if __name__ == "__main__":
    main()
