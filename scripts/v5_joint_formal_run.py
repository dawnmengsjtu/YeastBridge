#!/usr/bin/env python3
"""v5_joint_formal_run — V5 联合增益正式确认运行（预注册，单次）。

协议：docs/cross_species_match/V5_JOINT_FORMAL_PROTOCOL.md
配置：configs/v5_joint_formal.json（输入 SHA-256 冻结，运行前强制核对）
来源：match_diag_20260908/diag02_signal_audit.py A–U 段诊断（事后配置，已披露）。

用法：
  python scripts/v5_joint_formal_run.py \
    --config configs/v5_joint_formal.json --output <run_dir>
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cross_species_match import load_matrix


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def auc_p(pos, neg):
    u = mannwhitneyu(pos, neg, alternative="greater")
    return float(u.statistic / (len(pos) * len(neg))), float(u.pvalue)


def l2n(V):
    n = np.linalg.norm(V, axis=1, keepdims=True)
    return np.divide(V, n, out=np.zeros_like(V), where=n > 0)


def z_rows(R):
    Z = R - np.nanmean(R, 1, keepdims=True)
    sd = np.nanstd(R, 1, keepdims=True)
    Z = Z / np.where(sd > 0, sd, 1.0)
    return np.clip(np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0),
                   -50, 50).astype(np.float32)


def remove_pcs(V, npc):
    if npc == 0:
        return l2n(V)
    Vc = V - V.mean(axis=0, keepdims=True)
    _, _, Vt = np.linalg.svd(Vc, full_matrices=False)
    for i in range(npc):
        Vc = Vc - np.outer(Vc @ Vt[i], Vt[i])
    return l2n(Vc)


def signflip_p(d, seed, n=2048):
    obs = float(np.mean(d))
    hits = 0
    for t in range(n):
        s = np.random.default_rng(seed + t).choice([-1, 1], d.size)
        if float(np.mean(s * d)) >= obs:
            hits += 1
    return hits / n


def build_projection_maps(path):
    p = pd.read_csv(path, sep="\t", compression="gzip",
                    usecols=["species", "gene", "module_id", "weight"])
    mods_h = set(p[p.species == "human"].module_id)
    mods_y = set(p[p.species == "yeast"].module_id)
    shared = sorted(mods_h & mods_y)
    mod_idx = {m: i for i, m in enumerate(shared)}
    edges = {}
    for r in p.itertuples(index=False):
        i = mod_idx.get(r.module_id)
        if i is not None:
            edges.setdefault((r.species, r.gene), []).append((i, float(r.weight)))
    return edges, shared


def project_v2_matrix(X, col_genes, species, edges, n_modules):
    rows_, cols_, data_ = [], [], []
    for j, g in enumerate(col_genes):
        es = edges.get((species, g), ())
        if not es:
            continue
        tot = sum(abs(w) for _, w in es)
        for i, w in es:
            rows_.append(j)
            cols_.append(i)
            data_.append(w / tot)
    M = sp.csr_matrix((data_, (rows_, cols_)),
                      shape=(len(col_genes), n_modules))
    Xd = np.nan_to_num(np.asarray(X, dtype=np.float64), nan=0.0)
    num = Xd @ M
    den = np.abs(Xd) @ np.abs(M)
    V = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return l2n(V)


def es_space(term_gene_sets, idx_map, R):
    rows_, cols_ = [], []
    tlist = sorted(term_gene_sets)
    for j, t in enumerate(tlist):
        for g in term_gene_sets[t]:
            rows_.append(idx_map[g])
            cols_.append(j)
    B = sp.csr_matrix((np.ones(len(rows_), dtype=np.float32), (rows_, cols_)),
                      shape=(R.shape[1], len(tlist)))
    root = np.sqrt(np.array([len(term_gene_sets[t]) for t in tlist],
                            dtype=np.float32))
    return l2n((z_rows(R) @ B) / root[None, :]), root


def load_foundation_reprs(st, scf_npy, scyeast_prior_npy, raw_h5ad):
    """V5.2 双基座表示加载（2026-09-12 自 repr_grid 抽出，逻辑逐字不变）。

    人侧 scF：get_embedding.py（bulk, pre_normalized=F）逐行嵌入 K562 raw h5ad，
      按 bridge_setup 同一正则 keep 行对齐（断言与 normalized h5ad 行序一致）。
    酵母侧 scyeast：scyeast_embed.py --prior-space 的平坦输入先验基因空间 G(5812,200)；
      KO 签名 z_rows 投影：Fy = z_rows(yeast_matrix) @ B，B 按 systematic 对齐。
    返回 (Fh, Fy, cache_sha)，cache_sha 为三个冻结缓存文件的 SHA-256（写运行记录用）。
    """
    import re as _re
    import anndata as ad

    cfg = st["config"]
    root = Path(cfg["stage_root"])
    Ek1, Ey1 = st["Ek1"], st["Ey1"]

    # ---- 人侧 scF 表示 ----
    # keep 正则必须与 bridge_setup 逐字一致：(P1|P1P2)，不含纯 P2 行。
    # （2026-09-12 修复：原 repr_grid 误写 (P1|P1P2|P2)，10353 行 vs Ek1 9608 行，
    #   首次真跑即被行数断言拦下；bridge_setup 是冻结口径，以它为准。）
    pat = _re.compile(r"^\d+_([^_]+)_(P1|P1P2)_")
    a_raw = ad.read_h5ad(raw_h5ad)
    keep = np.array([bool(pat.match(str(n))) and
                     pat.match(str(n)).group(2) in {"P1", "P2", "P1P2"}
                     for n in a_raw.obs.index])
    Fh_all = np.load(scf_npy)
    assert Fh_all.shape[0] == a_raw.shape[0], \
        f"scF 嵌入行数 {Fh_all.shape[0]} != h5ad 行数 {a_raw.shape[0]}"
    a_norm = ad.read_h5ad(root / cfg["inputs"]["k562_h5ad"]["path"])
    assert list(a_raw.obs.index) == list(a_norm.obs.index), "raw 与 normalized h5ad 行序不一致"
    Fh = l2n(Fh_all[keep].astype(np.float64))
    assert Fh.shape[0] == Ek1.shape[0], f"Fh {Fh.shape} vs Ek1 {Ek1.shape} 行数不齐"

    # ---- 酵母侧 scyeast 先验基因空间表示 ----
    ym = pd.read_csv(root / cfg["inputs"]["yeast_matrix"]["path"], sep="\t", index_col=0)
    assert list(ym.index) == list(st["meta_y"].iloc[:, 0]) or len(ym) == Ey1.shape[0]
    G = np.load(scyeast_prior_npy).astype(np.float64)
    order = [l.strip() for l in open(root / "models/scyeast/gene_order_5812.txt")]
    pos = {g: i for i, g in enumerate(order)}
    B = np.zeros((ym.shape[1], G.shape[1]))
    for j, g in enumerate(ym.columns):
        if g in pos:
            B[j] = G[pos[g]]
    Fy = l2n(z_rows(ym.to_numpy()).astype(np.float64) @ B)
    assert Fy.shape[0] == Ey1.shape[0]

    cache_sha = {"scf_npy": sha256_file(scf_npy),
                 "scyeast_prior_npy": sha256_file(scyeast_prior_npy),
                 "raw_h5ad": sha256_file(raw_h5ad)}
    return Fh, Fy, cache_sha


def repr_grid(config_path, outdir, scf_npy, scyeast_prior_npy, raw_h5ad):
    """V5.2 表示层消融格（2026-09-11，事后变体，已与用户签）：
    {rawES, scF} × {rawES, scyeast 先验基因空间}，对齐器/拆分/row_negs/评测机器不变。
    双基座表示加载已抽出为 load_foundation_reprs（逐字同逻辑）。
    """
    st = bridge_setup(config_path)
    cfg = st["config"]
    Ek1, Ey1 = st["Ek1"], st["Ey1"]
    pairs = st["pairs"]
    ph, py_ = st["ph"], st["py"]
    pair_tr = np.array([p[2] == "train" for p in pairs])
    pair_te = np.array([p[2] == "test" for p in pairs])
    row_negs = st["row_negs"]
    ta = cfg["trackA"]

    Fh, Fy, _ = load_foundation_reprs(st, scf_npy, scyeast_prior_npy, raw_h5ad)

    def eval_cell(X, Y, label):
        Xtr, Ytr = X[ph[pair_tr]], Y[py_[pair_tr]]
        Mxy = Xtr.T @ Ytr
        best = (-1.0, None)
        for lam in ta["ridge_lambda_grid"]:
            W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(X.shape[1]), Mxy)
            S = (X @ W) @ Y.T
            a_tr, _ = auc_p(S[ph[pair_tr], py_[pair_tr]],
                            np.concatenate([S[r, row_negs[r]] for r in np.unique(ph[pair_tr])]))
            if a_tr > best[0]:
                best = (a_tr, W)
        W = best[1]
        S = (X @ W) @ Y.T
        auc, p = auc_p(S[ph[pair_te], py_[pair_te]],
                       np.concatenate([S[r, row_negs[r]] for r in np.unique(ph[pair_te])]))
        return {"repr": label, "test_auc": auc, "test_p": p,
                "ridge_train_auc": best[0], "dim_x": int(X.shape[1]), "dim_y": int(Y.shape[1])}

    cells = {
        "rawES_x_rawES(=V5.1口径)": (Ek1, Ey1),
        "scF_x_rawES": (Fh, Ey1),
        "rawES_x_scyeast_prior": (Ek1, Fy),
        "scF_x_scyeast_prior(双基座)": (Fh, Fy),
    }
    out = {"protocol": "V5.2 repr-grid (post-hoc variant, signed 2026-09-11)",
           "note": "对齐器/冻结拆分/row_negs/评测与 V5.1 同机; 注册参照: "
                   "armA=0.4948 armB=0.5169 armC=0.5687; ESM2/B2 蛋白路见 route_b 产物",
           "cells": {}}
    for label, (X, Y) in cells.items():
        r = eval_cell(X, Y, label)
        out["cells"][label] = r
        print(f"  {label}: test AUC={r['test_auc']:.4f} p={r['test_p']:.2g}", flush=True)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "repr_grid_result.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[done] -> {outdir/'repr_grid_result.json'}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--repr-grid", default=None,
                    help="V5.2 表示层消融格输出目录(需提供 scF/scyeast 缓存)")
    ap.add_argument("--scf-npy", default=None)
    ap.add_argument("--scyeast-prior-npy", default=None)
    ap.add_argument("--raw-h5ad", default=None)
    ap.add_argument("--export-gobridge-tasks", default=None,
                    help="逗号分隔人源靶点;导出 GO 线全库任务表(事后变体,不动正式协议)")
    ap.add_argument("--export-outdir", default=None)
    ap.add_argument("--export-repr", default="v51",
                    choices=["v51", "scf_scyeast"],
                    help="任务表导出的表示层:v51=冻结 rawES(注册口径,默认);"
                         "scf_scyeast=V5.2 双基座(事后变体,需 --scf-npy/"
                         "--scyeast-prior-npy/--raw-h5ad)")
    args = ap.parse_args()
    if args.export_gobridge_tasks:
        if not args.export_outdir:
            ap.error("--export-outdir is required with --export-gobridge-tasks")
        if args.export_repr == "scf_scyeast" and not (
                args.scf_npy and args.scyeast_prior_npy and args.raw_h5ad):
            ap.error("--export-repr scf_scyeast requires --scf-npy, "
                     "--scyeast-prior-npy and --raw-h5ad")
        export_gobridge_tasks(args.config,
                              [t.strip() for t in args.export_gobridge_tasks.split(",") if t.strip()],
                              args.export_outdir,
                              repr_variant=args.export_repr,
                              scf_npy=args.scf_npy,
                              scyeast_prior_npy=args.scyeast_prior_npy,
                              raw_h5ad=args.raw_h5ad)
        return
    if args.repr_grid:
        repr_grid(args.config, args.repr_grid, args.scf_npy,
                  args.scyeast_prior_npy, args.raw_h5ad)
        return
    if not args.output:
        ap.error("--output is required unless --repr-grid is given")
    cfg = json.loads(Path(args.config).read_text())
    root = Path(cfg["stage_root"])
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    res = {"protocol": cfg["protocol"], "config_sha256": sha256_file(args.config),
           "script_sha256": sha256_file(str(Path(__file__).resolve()))}

    # ---- 输入哈希核对（不匹配即中止）
    manifest = {}
    for k, v in cfg["inputs"].items():
        p = root / v["path"]
        h = sha256_file(str(p))
        manifest[k] = {"path": v["path"], "sha256": h, "match": h == v["sha256"]}
        if h != v["sha256"]:
            print(f"ABORT: input hash mismatch: {k} ({v['path']})")
            sys.exit(2)
    (out_dir / "INPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=1))
    print("input hashes: all match")

    import anndata as ad
    gate = cfg["data_quality_gate"]["column_abs_max_threshold"]

    # ---- 轨道 A 数据
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
    assert n_drop == cfg["data_quality_gate"]["expected_dropped_columns"], \
        f"column gate dropped {n_drop}, expected 73"
    Xp = Xp[:, col_ok]
    genes_k = [g for g, ok in zip(genes_k, col_ok) if ok]
    gidx_k = {g: j for j, g in enumerate(genes_k)}
    t2row = {}
    for j, t in enumerate(target_of):
        t2row.setdefault(t, j)
    print(f"K562 P1 profiles {Xp.shape[0]}, gate dropped {n_drop} columns")

    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    gidx_y = {g: j for j, g in enumerate(genes_y)}
    meta_y = pd.read_csv(root / cfg["inputs"]["yeast_metadata"]["path"], sep="\t")
    meta_y = meta_y[meta_y.benchmark_eligible.astype(str).str.lower() == "1"]
    meta_y = meta_y[meta_y.signature_id.isin(mat_y)].reset_index(drop=True)
    Y = np.stack([mat_y[s] for s in meta_y.signature_id])

    # goSlim 投影（轨道 B 用 + 轨道 A 的 valid_k 口径复现）
    edges_gs, shared_gs = build_projection_maps(
        root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)
    V_k = project_v2_matrix(Xp, genes_k, "human", edges_gs, N_MOD)
    valid_k = np.linalg.norm(V_k, axis=1) > 1e-9
    V_y_gs = project_v2_matrix(Y, genes_y, "yeast", edges_gs, N_MOD)

    # GO-BP ES 空间
    anc = pd.read_csv(root / cfg["inputs"]["go_bp_ancestors"]["path"],
                      sep="\t", compression="gzip",
                      usecols=["species", "gene", "ancestor_go_id"])
    tmin, tmax = cfg["representation"]["es_term_min_genes"], \
        cfg["representation"]["es_term_max_genes"]

    def term_sets(df):
        d = {}
        for g, t in zip(df.gene, df.ancestor_go_id):
            d.setdefault(t, []).append(g)
        return {t: gs for t, gs in d.items() if tmin <= len(gs) <= tmax}

    ts_h = term_sets(anc[(anc.species == "human") & anc.gene.isin(gidx_k)])
    ts_y = term_sets(anc[(anc.species == "yeast") & anc.gene.isin(gidx_y)])
    del anc
    terms = sorted(set(ts_h) & set(ts_y))
    print(f"GO-BP shared terms: {len(terms)}")
    E_k, root_h = es_space({t: ts_h[t] for t in terms}, gidx_k, Xp)
    E_y, _ = es_space({t: ts_y[t] for t in terms}, gidx_y, Y)
    E_k = E_k.astype(np.float64)
    E_y = E_y.astype(np.float64)
    Ek1 = remove_pcs(E_k, cfg["representation"]["pc_removal_for_trackA_and_members"])
    Ey1 = remove_pcs(E_y, cfg["representation"]["pc_removal_for_trackA_and_members"])

    # 同源对（官方 split）
    edges_df = pd.read_csv(root / cfg["inputs"]["orthodb_edges"]["path"], sep="\t")
    h_partner = {}
    for h, y in zip(edges_df.human_target_symbol, edges_df.yeast_deletion_target):
        h_partner.setdefault(h, set()).add(y)
    orf2row = {}
    for i, r in meta_y.iterrows():
        orf2row.setdefault(str(r.target_stable_id), i)
    pairs = []
    for _, e in edges_df.iterrows():
        h, y, spl = e.human_target_symbol, e.yeast_deletion_target, e.split
        j, k2 = t2row.get(h), orf2row.get(y)
        if j is None or k2 is None or not valid_k[j]:
            continue
        pairs.append((h, y, spl, j, k2))
    assert len(pairs) == 753, f"usable pairs {len(pairs)}, expected 753"
    ph = np.array([p[3] for p in pairs], dtype=int)
    py_ = np.array([p[4] for p in pairs], dtype=int)
    pair_tr = np.array([p[2] == "train" for p in pairs])
    pair_te = np.array([p[2] == "test" for p in pairs])
    ny = E_y.shape[0]
    ta = cfg["trackA"]

    # 查询行内负例（seed 冻结）
    sym_of_row = {}
    for p in pairs:
        sym_of_row.setdefault(p[3], p[0])
    row_negs = {}
    rng = np.random.default_rng(ta["seed_row_negatives"])
    for r in np.unique(ph):
        blk = {k2 for k2 in (orf2row.get(y) for y in h_partner.get(sym_of_row[r], ()))
               if k2 is not None}
        cand = np.setdiff1d(np.arange(ny), sorted(blk) or [-1])
        row_negs[r] = rng.choice(cand, size=min(ta["row_negatives_per_query"],
                                                len(cand)), replace=False)

    # 臂 C 的 ridge（λ 只按官方 train AUC 选）
    def fit_ridge():
        Xtr, Ytr = Ek1[ph[pair_tr]], Ey1[py_[pair_tr]]
        Mxy = Xtr.T @ Ytr
        best = (-1.0, None)
        for lam in ta["ridge_lambda_grid"]:
            W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
            S_w = (Ek1 @ W) @ Ey1.T
            a_tr, _ = auc_p(S_w[ph[pair_tr], py_[pair_tr]],
                            np.concatenate([S_w[r, row_negs[r]]
                                            for r in np.unique(ph[pair_tr])]))
            if a_tr > best[0]:
                best = (a_tr, W)
        return best[1], best[0]

    W_A, train_auc = fit_ridge()
    S_A = E_k @ E_y.T
    S_B = Ek1 @ Ey1.T

    def pls4_S():
        Xtr, Ytr = Ek1[ph[pair_tr]], Ey1[py_[pair_tr]]
        mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
        my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
        Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                     full_matrices=False)
        Px = l2n(((Ek1 - mx_) / sx_) @ Usv[:, :4])
        Py = l2n(((Ey1 - my_) / sy_) @ Vsv.T[:, :4])
        return Px @ Py.T

    def pct_rows(S):
        M = np.empty_like(S)
        for r in range(S.shape[0]):
            M[r] = (np.argsort(np.argsort(S[r])) + 1) / S.shape[1]
        return M

    recipe = ta.get("armC_recipe", {"type": "single_ridge"})
    resA = {}
    if recipe.get("type") == "rank_mean_ensemble":
        parts = []
        if "es_pc1_ridge" in recipe["members"]:
            parts.append((Ek1 @ W_A) @ Ey1.T)
        if "es_pc1_pls4" in recipe["members"]:
            parts.append(pls4_S())
        S_C = np.mean([pct_rows(S) for S in parts], axis=0)
    else:
        S_C = (Ek1 @ W_A) @ Ey1.T
    resA["armC_recipe"] = recipe
    resA["ridge_train_auc"] = train_auc
    te_rows = np.unique(ph[pair_te])

    def arm_eval(S, name):
        pos = S[ph[pair_te], py_[pair_te]]
        neg = np.concatenate([S[r, row_negs[r]] for r in te_rows])
        auc, p = auc_p(pos, neg)
        resA[f"{name}_test_auc"] = auc
        resA[f"{name}_test_p"] = p
        print(f"  A {name}: test auc={auc:.4f} p={p:.2g}")

    for nm, S in (("armA_no_align", S_A), ("armB_unsup_pc1", S_B),
                  ("armC_joint_ridge", S_C)):
        arm_eval(S, nm)

    # 配对 ΔAUC bootstrap（两臂共用同一重采样索引）
    def boot_delta(Sa, Sb, seed):
        pos_a = Sa[ph[pair_te], py_[pair_te]]
        neg_a = np.concatenate([Sa[r, row_negs[r]] for r in te_rows])
        pos_b = Sb[ph[pair_te], py_[pair_te]]
        neg_b = np.concatenate([Sb[r, row_negs[r]] for r in te_rows])
        wa = (pos_a[:, None] > neg_a[None, :]).mean(axis=1)
        wb = (pos_b[:, None] > neg_b[None, :]).mean(axis=1)
        rb = np.random.default_rng(seed)
        idx = np.arange(len(wa))
        ds = np.empty(ta["n_bootstrap"])
        for t in range(ta["n_bootstrap"]):
            s = rb.choice(idx, len(idx), replace=True)
            ds[t] = wa[s].mean() - wb[s].mean()
        return float(ds.mean()), float(np.percentile(ds, 2.5)), \
            float(np.percentile(ds, 97.5))

    sf_seed = ta["seed_signflip_base"]
    for nm, Sa, Sb in (("C_vs_A", S_C, S_A), ("C_vs_B", S_C, S_B)):
        m, lo, hi = boot_delta(Sa, Sb, ta["seed_bootstrap"])
        resA[f"boot_{nm}"] = {"mean": m, "ci95_low": lo, "ci95_high": hi}
        print(f"  A ΔAUC {nm}: {m:+.4f} CI[{lo:+.4f},{hi:+.4f}]")

    # OOF 5 折（人靶符号分组）
    sym_list = sorted(set(p[0] for p in pairs))
    fold_of_sym = {s: i % ta["oof_folds"] for i, s in enumerate(sym_list)}
    fold_of_pair = np.array([fold_of_sym[p[0]] for p in pairs])

    def pct_ranks(S, mask):
        return np.array([(np.sum(S[h_] > S[h_, y_]) + 1) / ny
                         for h_, y_ in zip(ph[mask], py_[mask])])

    oof = {nm: np.full(len(pairs), np.nan)
           for nm in ("A", "B", "C")}
    oof_hit = {nm: np.full(len(pairs), np.nan) for nm in ("B", "C")}
    for f in range(ta["oof_folds"]):
        m_tr, m_te = fold_of_pair != f, fold_of_pair == f
        Xtr, Ytr = Ek1[ph[m_tr]], Ey1[py_[m_tr]]
        Mxy = Xtr.T @ Ytr
        best = (-1.0, None)
        for lam in ta["ridge_lambda_grid"]:
            W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
            S_w = (Ek1 @ W) @ Ey1.T
            a_tr, _ = auc_p(S_w[ph[m_tr], py_[m_tr]],
                            np.concatenate([S_w[r, row_negs[r]]
                                            for r in np.unique(ph[m_tr])]))
            if a_tr > best[0]:
                best = (a_tr, W)
        Sc1 = (Ek1 @ best[1]) @ Ey1.T
        if recipe.get("type") == "rank_mean_ensemble":
            mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
            my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
            Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @
                                         ((Ytr - my_) / sy_), full_matrices=False)
            Sp = l2n(((Ek1 - mx_) / sx_) @ Usv[:, :4]) @ \
                l2n(((Ey1 - my_) / sy_) @ Vsv.T[:, :4]).T
            SC_f = np.mean([pct_rows(Sc1), pct_rows(Sp)], axis=0)
        else:
            SC_f = Sc1
        for nm, S in (("A", S_A), ("B", S_B), ("C", SC_f)):
            oof[nm][m_te] = pct_ranks(S, m_te)
        for i in np.flatnonzero(m_te):
            h_ = ph[i]
            partners = [orf2row[y2] for y2 in h_partner.get(sym_of_row[h_], ())
                        if orf2row.get(y2) is not None]
            if not partners:
                continue
            oof_hit["B"][i] = min(int(np.sum(S_B[h_] > S_B[h_, c_])) + 1
                                  for c_ in partners) <= 10
            oof_hit["C"][i] = min(int(np.sum(SC_f[h_] > SC_f[h_, c_])) + 1
                                  for c_ in partners) <= 10
    valid = ~np.isnan(oof["A"])
    for nm_a, nm_b in (("C", "A"), ("C", "B")):
        d = oof[nm_b][valid] - oof[nm_a][valid]
        sf_seed += 1
        p_ = signflip_p(d, sf_seed, ta["n_signflip"])
        resA[f"oof_{nm_a}_vs_{nm_b}"] = {"gain": float(d.mean()),
                                         "signflip_p": p_}
        print(f"  A OOF {nm_a}-vs-{nm_b}: {d.mean():+.4f} p={p_:.3g}")
    hb = oof_hit["B"][~np.isnan(oof_hit["B"])]
    hc = oof_hit["C"][~np.isnan(oof_hit["C"])]
    d_hit = hc.astype(float) - hb.astype(float)
    sf_seed += 1
    resA["any_hit10"] = {"rate_B": float(hb.mean()), "rate_C": float(hc.mean()),
                         "signflip_p": signflip_p(d_hit, sf_seed, ta["n_signflip"])}
    print(f"  A any-hit@10: B={hb.mean():.4f} C={hc.mean():.4f} "
          f"p={resA['any_hit10']['signflip_p']:.3g}")

    gateA_primary = (resA["armC_joint_ridge_test_p"] < 0.05 and
                     resA["boot_C_vs_A"]["ci95_low"] > 0)
    gateA_secondary = (resA["oof_C_vs_A"]["signflip_p"] < 0.05 or
                       resA["any_hit10"]["signflip_p"] < 0.05)
    resA["gate_primary_pass"] = bool(gateA_primary)
    resA["gate_secondary_at_least_one_pass"] = bool(gateA_secondary)
    print(f"  A gates: primary={gateA_primary} secondary={gateA_secondary}")

    pd.DataFrame({
        "human": [p[0] for p in pairs], "yeast_orf": [p[1] for p in pairs],
        "split": [p[2] for p in pairs], "fold": fold_of_pair,
        "oof_pct_rank_A": oof["A"], "oof_pct_rank_B": oof["B"],
        "oof_pct_rank_C": oof["C"]}).to_csv(
        out_dir / "per_query_ranks_trackA.tsv", sep="\t", index=False)
    res["trackA"] = resA

    # ---- 轨道 B 数据（Norman 基本结构）
    tb = cfg["trackB"]
    CLASSES = tb["label_universe"]
    CLS_IDX = {c: i for i, c in enumerate(CLASSES)}
    genes_n, mat_n = load_matrix(root / cfg["inputs"]["higher_matrix"]["path"])
    meta_n = pd.read_csv(root / cfg["inputs"]["higher_metadata"]["path"], sep="\t")
    meta_n = meta_n[meta_n.benchmark_eligible.astype(str).str.lower() == "true"]
    meta_n = meta_n[meta_n.signature_id.isin(mat_n)].reset_index(drop=True)
    assert len(meta_n) == 207, f"Norman sigs {len(meta_n)}, expected 207"
    Nn = np.stack([mat_n[s] for s in meta_n.signature_id])
    gidx_n = {g: j for j, g in enumerate(genes_n)}
    # Norman ES 完全复刻诊断配方：term 集 = 640 共享 term，基因取 K562 侧集合在
    # Norman 矩阵中存在者，归一用 K562 侧计数 root_h（与 diag02 B_n 构建一致）
    rows_n, cols_n = [], []
    for j_, t in enumerate(terms):
        for g in ts_h[t]:
            if g in gidx_n:
                rows_n.append(gidx_n[g])
                cols_n.append(j_)
    B_n = sp.csr_matrix((np.ones(len(rows_n), dtype=np.float32),
                         (rows_n, cols_n)), shape=(Nn.shape[1], len(terms)))
    E_n = l2n((z_rows(Nn) @ B_n) / root_h[None, :]).astype(np.float64)
    En1 = remove_pcs(E_n, 1)
    V_n = project_v2_matrix(Nn, genes_n, "human", edges_gs, N_MOD)

    lab_n = meta_n.function_label.values
    grp_n = meta_n.split_group.values
    te_n = np.array([("test" in str(s).lower()) for s in meta_n.split])
    in_cls = np.array([str(l) in CLS_IDX for l in lab_n])
    yeast_ref_mask = np.array([(str(m.split).lower() == "train" and
                                str(m.function_label) in CLS_IDX)
                               for _, m in meta_y.iterrows()])
    Yg0 = V_y_gs[yeast_ref_mask]
    Ey64r = E_y[yeast_ref_mask]
    Ey1r = Ey1[yeast_ref_mask]
    lab_y_ref = np.array([CLS_IDX[str(m.function_label)]
                          for _, m in meta_y[yeast_ref_mask].iterrows()])

    def cls_scores(qv, refs, ref_lab):
        sims = refs @ qv
        out = np.full(len(CLASSES), -np.inf)
        for c in range(len(CLASSES)):
            m = ref_lab == c
            if m.any():
                out[c] = sims[m].max()
        return out

    def rank_of_truth(scores, truth):
        if truth not in CLS_IDX:
            return 0
        order = np.argsort(-scores, kind="stable")
        return int(np.where(order == CLS_IDX[truth])[0][0]) + 1

    # 5 路成员分数
    S_h = np.full((len(meta_n), len(CLASSES)), -np.inf)
    for g in pd.unique(grp_n):
        href_m = (grp_n != g) & in_cls
        Href = En1[href_m]
        lab_h = np.array([CLS_IDX[str(x)] for x in lab_n[href_m]])
        for i in np.flatnonzero(grp_n == g):
            S_h[i] = cls_scores(En1[i], Href, lab_h)
    S_y_gs0 = np.stack([cls_scores(V_n[i], Yg0, lab_y_ref)
                        for i in range(len(meta_n))])
    S_jrepr = np.stack([cls_scores(En1[i] @ W_A, Ey1r, lab_y_ref)
                        for i in range(len(meta_n))])
    S_y_es0 = np.stack([cls_scores(E_n[i], Ey64r, lab_y_ref)
                        for i in range(len(meta_n))])

    def label_centroids(mat, lab_idx):
        return np.stack([mat[lab_idx == c].mean(0) if (lab_idx == c).any()
                         else np.zeros(mat.shape[1])
                         for c in range(len(CLASSES))])

    lam_lab = tb["label_centroid_ridge_lambda"]
    S_lab = np.full((len(meta_n), len(CLASSES)), -np.inf)
    for g in pd.unique(grp_n):
        trm = (grp_n != g) & in_cls
        h_lab = np.array([CLS_IDX[str(x)] for x in lab_n[trm]])
        Ch = label_centroids(V_n[trm], h_lab)
        Cy = label_centroids(Yg0, lab_y_ref)
        W = np.linalg.solve(Ch.T @ Ch + lam_lab * np.eye(Ch.shape[1]),
                            Ch.T @ Cy)
        for i in np.flatnonzero(grp_n == g):
            S_lab[i] = cls_scores(V_n[i] @ W, Yg0, lab_y_ref)

    MEMBERS = [S_h, S_y_gs0, S_jrepr, S_lab, S_y_es0]

    def rrf(k):
        out = np.zeros_like(MEMBERS[0], dtype=float)
        for S in MEMBERS:
            for i in range(S.shape[0]):
                out[i] += 1.0 / (k + np.argsort(np.argsort(-S[i])))
        return out

    def rr_of(Sf):
        rk = [rank_of_truth(Sf[i], str(lab_n[i])) for i in range(len(meta_n))]
        return np.array([1.0 / x if x > 0 else 0.0 for x in rk])

    rr_best_single = np.array([1.0 / x if x > 0 else 0.0
                               for x in [rank_of_truth(S_y_gs0[i], str(lab_n[i]))
                                         for i in range(len(meta_n))]])
    resB = {"single_arm_y_gs0_mrr": float(rr_best_single.mean())}
    # 嵌套选择（k 仅在其余组上选）
    rr_nest = np.full(len(meta_n), np.nan)
    picks = {}
    for g in pd.unique(grp_n):
        te = grp_n == g
        tr = (~te) & in_cls
        best_v, best_k = -1.0, None
        for k in tb["rrf_k_grid"]:
            Sf = rrf(k)
            v = np.mean([1.0 / x if x > 0 else 0.0 for x in
                         [rank_of_truth(Sf[i], str(lab_n[i]))
                          for i in np.flatnonzero(tr)]])
            if v > best_v:
                best_v, best_k = v, k
        picks[best_k] = picks.get(best_k, 0) + 1
        Sf = rrf(best_k)
        for i in np.flatnonzero(te):
            r_ = rank_of_truth(Sf[i], str(lab_n[i]))
            rr_nest[i] = 1.0 / r_ if r_ > 0 else 0.0
    sf_seed += 1
    d = rr_nest - rr_best_single
    resB["nested_joint_mrr"] = float(np.nanmean(rr_nest))
    resB["nested_delta_vs_best_single"] = float(np.nanmean(d))
    resB["nested_signflip_p"] = signflip_p(d[~np.isnan(d)], sf_seed,
                                           ta["n_signflip"])
    resB["nested_k_picks"] = {str(k): v for k, v in picks.items()}
    for k in (10, 20):
        rr_k = rr_of(rrf(k))
        resB[f"fixed_k{k}_mrr"] = float(rr_k.mean())
    resB["gate_pass"] = bool(resB["nested_signflip_p"] < 0.05)
    print(f"  B nested joint MRR={resB['nested_joint_mrr']:.3f} "
          f"Δ={resB['nested_delta_vs_best_single']:+.4f} "
          f"p={resB['nested_signflip_p']:.3g} gate={resB['gate_pass']}")

    pd.DataFrame({
        "signature_id": meta_n.signature_id, "label": lab_n,
        "split_group": grp_n, "official_test": te_n,
        "rr_best_single_y_gs0": rr_best_single, "rr_nested_joint": rr_nest,
    }).to_csv(out_dir / "per_query_rr_trackB.tsv", sep="\t", index=False)
    res["trackB"] = resB

    (out_dir / "result.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1))
    print("saved:", out_dir / "result.json")
    print(f"GATE A primary: {resA['gate_primary_pass']} | "
          f"GATE A secondary: {resA['gate_secondary_at_least_one_pass']} | "
          f"GATE B: {resB['gate_pass']}")


# ================================================================
# 可复用桥查询接口：供 yeastbridge_agent.py 等编排器调用。
# 复用 main() 的全部计算（ES 空间、对齐器、检索），不简化。
# ================================================================
import anndata as ad
from scipy.stats import mannwhitneyu


def _auc_p(pos, neg):
    u = mannwhitneyu(pos, neg, alternative="greater")
    return float(u.statistic / (len(pos) * len(neg))), float(u.pvalue)


def bridge_setup(config_path):
    """一次性加载冻结数据 → ES 空间 → 对齐器。返回可复用的桥状态。"""
    cfg = json.loads(Path(config_path).read_text())
    root = Path(cfg["stage_root"])
    gate = cfg["data_quality_gate"]["column_abs_max_threshold"]

    a = ad.read_h5ad(root / cfg["inputs"]["k562_h5ad"]["path"])
    genes_k = list(a.var.gene_name)
    obs_names = list(a.obs.index)
    pat = re.compile(r"^\d+_([^_]+)_(P1|P1P2)_")
    keep, target_of = [], []
    for name in obs_names:
        m = pat.match(name)
        if m and m.group(2) in {"P1", "P2", "P1P2"}:
            keep.append(name); target_of.append(m.group(1))
    Xp = a.X[[obs_names.index(n) for n in keep]]
    Xp = Xp.toarray() if hasattr(Xp, "toarray") else np.asarray(Xp)
    Xp = np.nan_to_num(Xp, nan=0.0)
    col_ok = np.abs(Xp).max(0) < gate
    Xp = Xp[:, col_ok]
    genes_k = [g for g, ok in zip(genes_k, col_ok) if ok]
    gidx_k = {g: j for j, g in enumerate(genes_k)}
    t2row = {}
    for j, t in enumerate(target_of):
        t2row.setdefault(t, j)

    genes_y, mat_y = load_matrix(root / cfg["inputs"]["yeast_matrix"]["path"])
    gidx_y = {g: j for j, g in enumerate(genes_y)}
    meta_y = pd.read_csv(root / cfg["inputs"]["yeast_metadata"]["path"], sep="\t")
    meta_y = meta_y[meta_y.benchmark_eligible.astype(str).str.lower() == "1"]
    meta_y = meta_y[meta_y.signature_id.isin(mat_y)].reset_index(drop=True)
    Y = np.stack([mat_y[s] for s in meta_y.signature_id])
    orf2row = {}
    for i, r in meta_y.iterrows():
        orf2row.setdefault(str(r.target_stable_id), i)

    edges, shared_gs = build_projection_maps(
        root / cfg["inputs"]["goslim_projection"]["path"])
    N_MOD = len(shared_gs)

    anc = pd.read_csv(root / cfg["inputs"]["go_bp_ancestors"]["path"],
                      sep="\t", compression="gzip",
                      usecols=["species", "gene", "ancestor_go_id"])
    tmin = cfg["representation"]["es_term_min_genes"]
    tmax = cfg["representation"]["es_term_max_genes"]

    def _term_sets(df):
        d = {}
        for g, t in zip(df.gene, df.ancestor_go_id):
            d.setdefault(t, []).append(g)
        return {t: gs for t, gs in d.items() if tmin <= len(gs) <= tmax}

    ts_h = _term_sets(anc[(anc.species == "human") & anc.gene.isin(gidx_k)])
    ts_y = _term_sets(anc[(anc.species == "yeast") & anc.gene.isin(gidx_y)])
    del anc
    terms = sorted(set(ts_h) & set(ts_y))

    def _es(ts, idx, R):
        rows_, cols_ = [], []
        for j, t in enumerate(terms):
            for g in ts[t]:
                rows_.append(idx[g]); cols_.append(j)
        B = sp.csr_matrix((np.ones(len(rows_), dtype=np.float32), (rows_, cols_)),
                          shape=(R.shape[1], len(terms)))
        root_arr = np.sqrt(np.array([len(ts[t]) for t in terms], dtype=np.float32))
        return l2n((z_rows(R) @ B) / root_arr[None, :]), root_arr

    E_k, _ = _es({t: ts_h[t] for t in terms}, gidx_k, Xp)
    E_y, _ = _es({t: ts_y[t] for t in terms}, gidx_y, Y)
    E_k = E_k.astype(np.float64)
    E_y = E_y.astype(np.float64)
    Ek1 = remove_pcs(E_k, cfg["representation"]["pc_removal_for_trackA_and_members"])
    Ey1 = remove_pcs(E_y, cfg["representation"]["pc_removal_for_trackA_and_members"])

    edges_df = pd.read_csv(root / cfg["inputs"]["orthodb_edges"]["path"], sep="\t")
    h_partner = {}
    for h, y in zip(edges_df.human_target_symbol, edges_df.yeast_deletion_target):
        h_partner.setdefault(h, set()).add(y)
    pairs = []
    for _, e in edges_df.iterrows():
        j, k2 = t2row.get(e.human_target_symbol), orf2row.get(e.yeast_deletion_target)
        if j is not None and k2 is not None:
            pairs.append((j, k2, e.split))
    ph = np.array([p[0] for p in pairs])
    py_ = np.array([p[1] for p in pairs])
    pair_tr = np.array([p[2] == "train" for p in pairs])

    row_negs = {}
    rng = np.random.default_rng(cfg["trackA"]["seed_row_negatives"])
    sym_of_row = {j: target_of[j] for j in ph}
    for r in np.unique(ph):
        blk = {k2 for k2 in (orf2row.get(y) for y in h_partner.get(sym_of_row[r], ()))
               if k2 is not None}
        cand = np.setdiff1d(np.arange(len(meta_y)), sorted(blk) or [-1])
        row_negs[r] = rng.choice(cand, size=min(cfg["trackA"]["row_negatives_per_query"],
                                                len(cand)), replace=False)

    # 对齐器：ridge λ 网格（train AUC 选参）+ PLS4
    Xtr, Ytr = Ek1[ph[pair_tr]], Ey1[py_[pair_tr]]
    Mxy = Xtr.T @ Ytr
    best = (-1.0, None)
    for lam in cfg["trackA"]["ridge_lambda_grid"]:
        W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
        S = (Ek1 @ W) @ Ey1.T
        auc, _ = _auc_p(S[ph[pair_tr], py_[pair_tr]],
                        np.concatenate([S[r, row_negs[r]] for r in np.unique(ph[pair_tr])]))
        if auc > best[0]:
            best = (auc, W)
    W_ridge = best[1]
    mx_ = Xtr.mean(0); sx_ = np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_ = Ytr.mean(0); sy_ = np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                 full_matrices=False)
    Wx4 = Usv[:, :4]; Wy4 = Vsv.T[:, :4]

    return {
        "config": cfg, "Ek1": Ek1, "Ey1": Ey1, "W_ridge": W_ridge,
        "Wx4": Wx4, "Wy4": Wy4, "mx": mx_, "sx": sx_, "my": my_, "sy": sy_,
        "orf2row": orf2row, "t2row": t2row, "h_partner": h_partner,
        "pairs": pairs, "ph": ph, "py": py_, "pair_tr": pair_tr,
        "row_negs": row_negs, "genes_y": genes_y, "meta_y": meta_y,
        "ridge_train_auc": best[0], "n_terms": len(terms),
    }


def bridge_query(target, state):
    """对单个靶点运行 V5.1 桥：投影 → 对齐 → 检索 → 返回排名。"""
    j = state["t2row"].get(target)
    if j is None:
        return {"status": "NO_GO", "reason": "靶点无 K562 CRISPRi 签名"}
    Ey1 = state["Ey1"]; Ek1 = state["Ek1"]
    W = state["W_ridge"]
    r1 = Ey1 @ (Ek1[j] @ W)
    p4 = ((Ek1[j] - state["mx"]) / state["sx"]) @ state["Wx4"]
    p4 = p4 / (np.linalg.norm(p4) if np.linalg.norm(p4) > 0 else 1)
    r2 = l2n(((Ey1 - state["my"]) / state["sy"]) @ state["Wy4"]) @ p4
    pct = lambda v: (np.argsort(np.argsort(v)) + 1) / len(v)
    ens = (pct(r1) + pct(r2)) / 2.0  # V5.1 C 臂集成（ridge+PLS4 秩平均）

    partners = state["h_partner"].get(target, set())
    results = []
    for y in partners:
        k2 = state["orf2row"].get(y)
        if k2 is None:
            continue
        rank = int(np.sum(ens > ens[k2])) + 1
        results.append({"ortholog": y, "rank": rank,
                        "pool": len(ens), "pct": rank / len(ens),
                        "ensemble_score": float(ens[k2])})
    return {"status": "PASS", "results": results,
            "method": "V5.1 C-arm (ridge+PLS4 rank ensemble)"}


def export_gobridge_tasks(config_path, targets, outdir,
                          repr_variant="v51", scf_npy=None,
                          scyeast_prior_npy=None, raw_h5ad=None,
                          state=None):
    """V5.2 GO 线任务表导出（2026-09-12 事后变体）。

    state: 可传入已加载的 bridge_setup 状态（总脚本 Phase 2 / 批量模式复用，
    避免逐靶点重复加载）；None 时自行 bridge_setup，行为与原版逐字不变。

    对指定人源靶点，用与 bridge_query 逐字相同的 V5.1 C 臂集成分
    （ridge + PLS4 秩平均）计算全库酵母基因排名，导出 yeast_task_<T>.tsv
    （rank, yeast_gene, ensemble_score），供 product_execute_hiphop.py 作为
    GO 线 task axis 与 B2 线同台对比。与 bridge_query 的唯一差异是输出：
    不按 h_partner（同源伙伴）截断，返回全库排名。不触碰正式协议路径。

    repr_variant（2026-09-12 增，GO 线效应量修复杠杆）：
      "v51"（默认）＝冻结 rawES 表示，注册口径，行为逐字不变；
      "scf_scyeast" ＝表示层换成 V5.2 双基座（人侧 scF × 酵母侧 scyeast
      先验空间，load_foundation_reprs），对齐配方仍为与 bridge_setup 逐字
      同构的 ridge λ 网格（train AUC 选参）+ PLS4 秩平均，在冻结 train
      pairs 上重训（与 repr_grid 消融格同配方），不触碰 test pairs。
      缓存文件 SHA-256 与重训 ridge_train_auc 写入 export_record.json。
    """
    if repr_variant not in ("v51", "scf_scyeast"):
        raise ValueError(f"unknown repr_variant: {repr_variant}")
    st = state if state is not None else bridge_setup(config_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    orfs = [str(x) for x in st["meta_y"]["target_stable_id"]]
    pct = lambda v: (np.argsort(np.argsort(v)) + 1) / len(v)
    report = {"mode": "export_gobridge_tasks",
              "repr_variant": repr_variant,
              "config": str(config_path),
              "config_sha256": sha256_file(str(config_path)),
              "script_sha256": sha256_file(str(Path(__file__).resolve())),
              "targets": {}}
    if repr_variant == "scf_scyeast":
        X, Y, cache_sha = load_foundation_reprs(
            st, scf_npy, scyeast_prior_npy, raw_h5ad)
        report["foundation_cache_sha256"] = cache_sha
        # ---- 对齐器重训：与 bridge_setup 逐字同配方，仅表示不同 ----
        cfg = st["config"]
        ph, py_ = st["ph"], st["py"]
        pair_tr = st["pair_tr"]
        row_negs = st["row_negs"]
        Xtr, Ytr = X[ph[pair_tr]], Y[py_[pair_tr]]
        Mxy = Xtr.T @ Ytr
        best = (-1.0, None)
        for lam in cfg["trackA"]["ridge_lambda_grid"]:
            W = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Mxy)
            S = (X @ W) @ Y.T
            a_tr, _ = auc_p(S[ph[pair_tr], py_[pair_tr]],
                            np.concatenate([S[r, row_negs[r]]
                                            for r in np.unique(ph[pair_tr])]))
            if a_tr > best[0]:
                best = (a_tr, W)
        W = best[1]
        mx_ = Xtr.mean(0); sx_ = np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
        my_ = Ytr.mean(0); sy_ = np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
        Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                    full_matrices=False)
        Wx4 = Usv[:, :4]; Wy4 = Vsv.T[:, :4]
        report["ridge_train_auc"] = best[0]
        print(f"  [scf_scyeast] refit aligner: ridge_train_auc={best[0]:.4f}",
              flush=True)
    else:
        X = st["Ek1"]; Y = st["Ey1"]; W = st["W_ridge"]
        mx_, sx_, my_, sy_ = st["mx"], st["sx"], st["my"], st["sy"]
        Wx4, Wy4 = st["Wx4"], st["Wy4"]
    report["pool"] = int(Y.shape[0])
    for t in targets:
        j = st["t2row"].get(t)
        if j is None:
            report["targets"][t] = {"status": "NO_K562_SIGNATURE"}
            print(f"  {t}: NO_K562_SIGNATURE", flush=True)
            continue
        r1 = Y @ (X[j] @ W)
        p4 = ((X[j] - mx_) / sx_) @ Wx4
        p4 = p4 / (np.linalg.norm(p4) if np.linalg.norm(p4) > 0 else 1)
        r2 = l2n(((Y - my_) / sy_) @ Wy4) @ p4
        ens = (pct(r1) + pct(r2)) / 2.0  # 与 bridge_query 相同的 C 臂集成
        order = np.argsort(-ens)
        seen, rows = set(), []
        for r in order:
            orf = orfs[int(r)]
            if orf in seen:
                continue
            seen.add(orf)
            rows.append((len(rows) + 1, orf, float(ens[int(r)])))
        fp = outdir / f"yeast_task_{t}.tsv"
        with fp.open("w", encoding="utf-8", newline="") as fh:
            fh.write("rank\tyeast_gene\tensemble_score\n")
            for rank, orf, sc in rows:
                fh.write(f"{rank}\t{orf}\t{sc:.6f}\n")
        report["targets"][t] = {"status": "OK", "rows": len(rows), "path": str(fp)}
        print(f"  {t}: {len(rows)} 行 -> {fp}", flush=True)
    (outdir / "export_record.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"[done] -> {outdir / 'export_record.json'}", flush=True)


if __name__ == "__main__":
    main()
