#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_ablations.py — 消融实验(plan §七 高风险项 + T3 混杂校正)

消融1(路线A 无同源基因退化, plan 高风险):
    T2/T3 按 routeA_init.tsv 的 init_type 分组(ortholog-init 2431 / random-init 4302)。
    假设:路线A 在 random-init 基因上退化;路线B(蛋白桥)天然免疫,应无此落差。
消融2(T3 混杂/批次校正, plan 高风险):
    (a) mean-profile 残差化:预测与真值同减训练集均值谱后重算 Spearman——检验"超出共享成分"的真实信号;
    (b) 跨子系列泛化:GSE42526 train→GSE42527 test 及反向——批次迁移;
    (c) responsive_class 子集:responsive 突变体(700)vs non_responsive(784)分组 Spearman——
        强响应子集才是模型可展示信号的区间(文献 O'Duibhir 2014 生长速率混杂的同族处理)。

用法:
    cd /public/home/mengxl/dzy/yeastbridge
    /public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/routeA/run_ablations.py

输出:
    results/routeA/ablation_no_ortholog.json
    results/harness/ablation_t3_confound.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[3]
ASSETS = Path("/public/home/mengxl/dzy/yeastbridge")  # legacy frozen asset root (scGPT weights/h5ad/mappings/external src); outputs go to this project
sys.path.insert(0, str(ROOT))
from eval import data  # noqa: E402

SEED = 42
FEATURES = ["routea_scgpt_ft", "routeb_protein", "esm2_mean", "routed_scyeast"]
ALPHAS = (0.1, 1.0, 10.0, 100.0)


def load_init_map():
    df = pd.read_csv(ROOT / "results_model_selection" / "step4_route_confirmation" / "routeA_assets" / "routeA_init.tsv", sep="\t", dtype=str)
    return dict(zip(df["systematic"], df["init_type"]))


def fit_ridge_auto(X, Y, tr, va):
    best = None
    for a in ALPHAS:
        m = Ridge(alpha=a).fit(X[tr], Y[tr])
        Pv = m.predict(X[va])
        sv = np.nanmean([spearmanr(Pv[i], Y[va][i]).statistic for i in range(len(Pv))])
        if best is None or sv > best[1]:
            best = (a, sv, m)
    return best[0], best[2]


def spear_rows(P, T):
    return np.array([spearmanr(P[i], T[i]).statistic for i in range(len(P))])


def t3_aligned(feature):
    split, expr = data.kemmeren()
    mut_sys = split.set_index("mutant_name")["systematic"].to_dict()
    emb_genes = data.feature_genes(feature)
    rows = []
    for m in expr["mutant"]:
        g = mut_sys.get(m, "")
        rows.append(g if g in emb_genes else None)
    expr = expr.assign(_ko=rows)
    expr = expr[expr["_ko"].notna()].reset_index(drop=True)
    X, kept = data.feature_matrix(feature, expr["_ko"].tolist())
    expr = expr[expr["_ko"].isin(kept)].reset_index(drop=True)
    gene_cols = [c for c in expr.columns if c not in ("mutant", "_ko")]
    Y = expr[gene_cols].to_numpy(np.float32)
    meta = split.set_index("mutant_name")
    return (X, Y, expr,
            expr["mutant"].map(meta["split"]).to_numpy(),
            expr["mutant"].map(meta["subseries"]).to_numpy(),
            expr["mutant"].map(meta["responsive_class"]).to_numpy())


def ablation_no_ortholog(init_map):
    out = {"t2": {}, "t3": {}, "group_def": "routeA_init.tsv init_type(ortholog=2431, random=4302)"}
    # --- T2 分组 ---
    lab = data.essentiality()
    lab = lab[lab["essentiality"].isin(["essential", "nonessential"])].reset_index(drop=True)
    for feat in FEATURES:
        X, kept = data.feature_matrix(feat, lab["systematic"].tolist())
        lab_f = lab[lab["systematic"].isin(kept)].reset_index(drop=True)
        y = (lab_f["essentiality"] == "essential").astype(int).to_numpy()
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, stratify=y, random_state=SEED)
        clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(Xtr, ytr)
        p = clf.predict_proba(Xte)[:, 1]
        genes_te = lab_f["systematic"].to_numpy()[
            train_test_split(np.arange(len(y)), test_size=0.2, stratify=y, random_state=SEED)[1]]
        row = {}
        for grp in ("ortholog", "random"):
            m = np.array([init_map.get(g) == grp for g in genes_te])
            row[grp] = {"auroc": float(roc_auc_score(yte[m], p[m])),
                        "auprc": float(average_precision_score(yte[m], p[m])),
                        "n_test": int(m.sum())}
        out["t2"][feat] = row
        print(f"[消融1-T2] {feat}: ortholog AUROC {row['ortholog']['auroc']:.3f} "
              f"vs random {row['random']['auroc']:.3f}")
    # --- T3 分组(按被敲基因 init_type) ---
    for feat in FEATURES:
        X, Y, expr, sl, sub, rcl = t3_aligned(feat)
        tr, va, te = sl == "train", sl == "val", sl == "test"
        a, mdl = fit_ridge_auto(X, Y, tr, va)
        sp = spear_rows(mdl.predict(X[te]), Y[te])
        grp_of = np.array([init_map.get(g, "random") for g in expr["_ko"].to_numpy()])[te]
        row = {}
        for grp in ("ortholog", "random"):
            m = grp_of == grp
            row[grp] = {"spearman_mean": float(np.nanmean(sp[m])), "n_test": int(m.sum())}
        row["alpha"] = a
        out["t3"][feat] = row
        print(f"[消融1-T3] {feat}: ortholog {row['ortholog']['spearman_mean']:.3f} "
              f"vs random {row['random']['spearman_mean']:.3f} (α={a})")
    return out


def ablation_confound():
    out = {"residualized": {}, "cross_subseries": {}, "responsive_subset": {}}
    for feat in FEATURES:
        X, Y, expr, sl, sub, rcl = t3_aligned(feat)
        tr, va, te = sl == "train", sl == "val", sl == "test"
        a, mdl = fit_ridge_auto(X, Y, tr, va)
        P, T = mdl.predict(X[te]), Y[te]
        mu = Y[tr].mean(0, keepdims=True)
        sp_raw = spear_rows(P, T)
        sp_res = spear_rows(P - mu, T - mu)
        out["residualized"][feat] = {"raw": float(np.nanmean(sp_raw)),
                                     "residual": float(np.nanmean(sp_res)), "alpha": a}
        print(f"[消融2-残差化] {feat}: raw {np.nanmean(sp_raw):.3f} -> residual {np.nanmean(sp_res):.3f}")
        # 跨子系列(α 沿用本特征 val 选出的)
        cs = {}
        for s_tr, s_te in (("GSE42526", "GSE42527"), ("GSE42527", "GSE42526")):
            m2 = Ridge(alpha=a).fit(X[(sl == "train") & (sub == s_tr)], Y[(sl == "train") & (sub == s_tr)])
            mte = (sub == s_te) & te
            cs[f"{s_tr}->{s_te}"] = float(np.nanmean(spear_rows(m2.predict(X[mte]), Y[mte])))
        out["cross_subseries"][feat] = cs
        print(f"[消融2-跨批次] {feat}: {cs}")
        rcl_te = rcl[te]
        rs = {}
        for cls in ("responsive", "non_responsive"):
            rs[cls] = {"spearman_mean": float(np.nanmean(sp_raw[rcl_te == cls])),
                       "n_test": int((rcl_te == cls).sum())}
        out["responsive_subset"][feat] = rs
        print(f"[消融2-响应子集] {feat}: responsive {rs['responsive']['spearman_mean']:.3f} "
              f"vs non_responsive {rs['non_responsive']['spearman_mean']:.3f}")
    return out


def main():
    init_map = load_init_map()
    a1 = ablation_no_ortholog(init_map)
    p1 = ROOT / "results_model_selection" / "step4_route_confirmation" / "ablations" / "ablation_no_ortholog.json"
    p1.write_text(json.dumps(a1, indent=1, ensure_ascii=False))
    a2 = ablation_confound()
    p2 = ROOT / "results_model_selection" / "step4_route_confirmation" / "ablations" / "ablation_t3_confound.json"
    p2.write_text(json.dumps(a2, indent=1, ensure_ascii=False))
    print(f"\n[out] {p1}\n[out] {p2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
