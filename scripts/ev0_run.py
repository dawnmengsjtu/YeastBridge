#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EV0 — 检索层外样本（第二采样确认）。协议：EV_PROTOCOL.md §EV0（先冻结后运行）。

对集：OrthoDB(ORF 行) ∪ OMA ∪ InParanoid − 753 官方对（对级排除）；
查询限于有 ESM-2 嵌入的人基因（810 表）、酵母限于池内（交集运行，如实登记）。
算子：V8 冻结配方——臂 A 裸余弦；臂 B 去 PC1（fx/fy 在冻结资产上拟合）；
臂 C 用冻结超参（npc*=3, ridge=1.0, anchored=1.0, pls=8，断言与
esm2_joint_tasks.json 一致）在官方 753 train 对上确定性重拟合后原样应用。
判据（冻结）：主门 = 臂 C MWU p<0.05 且 ΔAUC(C−A) bootstrap CI95 下界>0。
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path("/public/home/mengxl/dzy/YeastBridge_0912/scripts")
sys.path.insert(0, str(HERE))
from v5_joint_formal_run import auc_p, l2n, sha256_file  # noqa: E402
from v7_b2_formal_run import pct_rows, v7_setup  # noqa: E402
YB_ROOT = Path("/public/home/mengxl/dzy/YeastBridge_0912")
from v8_esm2_formal_run import apply_pc, fit_pc  # noqa: E402

MAP = Path("/public/home/mengxl/dzy/YeastBridge_0912/raw/externalvalidation/mappings")
EV = Path("/public/home/mengxl/dzy/YeastBridge_0912")
IS_X = (EV / "raw/externalvalidation/ev0x_embed_merged/esm2_mean_fp32.npy").exists()
OUT = EV / ("results_a6000/externalvalidation/ev0x_retrieval" if IS_X else "results_a6000/externalvalidation/ev0_retrieval")
OUT.mkdir(parents=True, exist_ok=True)
ORF = re.compile(r"^Y[A-Z]\d{2}[0-9A-Z]?[A-Z]?$")
NPC_DIRECT, NPC_ALIGNED, LAM_R, LAM_A, PLS_K = 1, 3, 1.0, 1.0, 8


def main():
    res = {"mode": "ev0_retrieval", "protocol": "EV0"}
    # ---- 冻结超参断言（与 EJ 方法冻结配置一致）----
    ej = json.loads((HERE.parent / "configs_a6000_frozen/esm2_joint_tasks.json")
                    .read_text())["method"]
    assert (ej["npc_direct"], ej["npc_aligned"], ej["ridge_lambda"],
            ej["anchored_lambda"], ej["pls_k"]) == (
        NPC_DIRECT, NPC_ALIGNED, LAM_R, LAM_A, PLS_K), "EV0 超参与冻结配置不符"
    res["hyper_frozen"] = dict(npc_direct=NPC_DIRECT, npc_aligned=NPC_ALIGNED,
                               ridge=LAM_R, anchored=LAM_A, pls_k=PLS_K)

    st = v7_setup(str(HERE.parent / "configs_a6000_frozen/v7_b2_joint_formal.json"))
    sym2row, pool_pos = st["sym2row_h"], st["pool_pos"]
    Y = st["Y_es"]
    Xall = st["X_es_all"]
    fx, fy = fit_pc(Xall), fit_pc(Y)
    prim753 = st["prim"]
    res["n_pairs_753_available"] = len(prim753)
    p753 = {(h, y) for h, y, _, _ in prim753}
    train_genes_753 = {h for h, _, s, _ in prim753}  # 全部 753 基因（含 test 折）

    # ---- 三源并集 ----
    pairs_src = {}
    odb = pd.read_csv(MAP / "orthodb_yeast_human.tsv", sep="\t")
    n_odb_raw = len(odb)
    odb = odb[odb.yeast_gene.astype(str).str.match(ORF)]
    for y, h in zip(odb.yeast_gene, odb.human_gene):
        pairs_src.setdefault((str(h).upper(), str(y)), set()).add("orthodb")
    oma = pd.read_csv(MAP / "oma_yeast_human.tsv", sep="\t")
    for y, h in zip(oma.yeast_systematic, oma.human_symbol):
        pairs_src.setdefault((str(h).upper(), str(y)), set()).add("oma")
    inp = pd.read_csv(MAP / "inparanoid_yeast_human.tsv", sep="\t")
    for y, h in zip(inp.yeast_systematic, inp.human_symbol):
        pairs_src.setdefault((str(h).upper(), str(y)), set()).add("inparanoid")
    res["union_pairs_raw"] = {"orthodb_rows": int(n_odb_raw),
                              "orthodb_orf_kept": len(odb),
                              "union_total": len(pairs_src)}

    # ---- 查询嵌入源扩展（EV0x）：优先用补推理合并表，回退 universe ----
    merged = EV / "raw/externalvalidation/ev0x_embed_merged"
    res["inputs_sha"] = {}
    if (merged / "esm2_mean_fp32.npy").exists():
        U = np.load(merged / "esm2_mean_fp32.npy").astype(np.float64)
        uidx = pd.read_csv(merged / "index.tsv", sep="\t").fillna("")
        res["inputs_sha"]["ev0x_embed"] = sha256_file(str(merged / "esm2_mean_fp32.npy"))
    else:
        uni = sorted((YB_ROOT / "raw/tier1_esm2").glob("universe*/"))
        udir = uni[0]
        U = np.load(udir / "esm2_mean_fp32.npy").astype(np.float64)
        uidx = pd.read_csv(udir / "index.tsv", sep="	").fillna("")
        res["inputs_sha"]["universe_esm2"] = sha256_file(str(udir / "esm2_mean_fp32.npy"))
    usyms = [str(s) for s in uidx["common"]]
    extra = [i for i, s in enumerate(usyms)
             if s and s not in sym2row and s.upper() not in sym2row]
    Xext = l2n(U[extra])
    ext_row = {usyms[i].upper(): len(sym2row) + k for k, i in enumerate(extra)}
    sym2row = {**{k.upper(): v for k, v in sym2row.items()}, **ext_row}
    Xall_q = np.vstack([Xall, Xext])
    res["query_embed_source"] = {"table810": int(Xall.shape[0]),
                                 "extended_added": len(extra),
                                 "source": "ev0x_merged" if (merged / "esm2_mean_fp32.npy").exists() else "universe"}

    cand = [(h, y) for (h, y) in pairs_src
            if h in sym2row and y in pool_pos and (h, y) not in p753]
    # 查询基因的全部并集伙伴（负例剔除用）
    partners = {}
    for h, y in pairs_src:
        if y in pool_pos:
            partners.setdefault(h, set()).add(y)
    qgenes = sorted({h for h, _ in cand})
    res["ev0_pairs"] = len(cand)
    res["ev0_query_genes"] = len(qgenes)
    res["gene_level_sensitivity"] = {
        "queries_seen_in_753": sorted(g for g in qgenes if g in train_genes_753)}

    Xq = Xall_q[[sym2row[g] for g in qgenes]]
    Y1 = apply_pc(Y, fy[0], fy[1], NPC_DIRECT)
    XN = apply_pc(Xq, fx[0], fx[1], NPC_ALIGNED)
    YN = apply_pc(Y, fy[0], fy[1], NPC_ALIGNED)

    # ---- 臂 C 冻结算子：官方 train 对重拟合（确定性，与 V8 逐位等价）----
    q_pos = st["q_pos"]
    qr = np.array([q_pos[h] for h, _, _, _ in prim753])
    yr = np.array([pool_pos[y] for _, y, _, _ in prim753])
    tr = np.array([s == "train" for _, _, s, _ in prim753])
    hrow753 = [sym2row[g] for g in st["query_genes"]]
    X753 = apply_pc(Xall[hrow753], fx[0], fx[1], NPC_ALIGNED)
    Xtr, Ytr = X753[qr[tr]], YN[yr[tr]]
    XtX, XtY = Xtr.T @ Xtr, Xtr.T @ Ytr
    I = np.eye(Xtr.shape[1])
    W_r = np.linalg.solve(XtX + LAM_R * I, XtY)
    W_a = np.linalg.solve(XtX + LAM_A * I, XtY + LAM_A * I)
    mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    U, _, Vt = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                             full_matrices=False)
    Px = l2n(((XN - mx_) / sx_) @ U[:, :PLS_K])
    Py = l2n(((YN - my_) / sy_) @ Vt[:PLS_K].T)
    S_r = (XN @ W_r) @ YN.T
    S_a = (XN @ W_a) @ YN.T
    S_p = Px @ Py.T
    S_Cq = np.mean([pct_rows(S_r), pct_rows(S_a), pct_rows(S_p)], axis=0)

    S_A = Xq @ Y.T
    S_B = apply_pc(Xq, fx[0], fx[1], 1) @ Y1.T

    # ---- 冻结负例（seed=20260916，排除该查询全部并集伙伴）----
    rng = np.random.default_rng(20260916)
    negs = {}
    n_pool = Y.shape[0]
    for g in qgenes:
        blk = sorted(partners.get(g, set()))
        candn = np.setdiff1d(np.arange(n_pool), blk or [-1])
        negs[g] = rng.choice(candn, size=min(30, len(candn)), replace=False)
    qi = {g: i for i, g in enumerate(qgenes)}
    pos_idx = [(qi[h], pool_pos[y]) for h, y in cand]

    def ev(S, name):
        pos = np.array([S[a, b] for a, b in pos_idx])
        neg = np.concatenate([S[qi[g], negs[g]] for g in qgenes])
        a, p = auc_p(pos, neg)
        res[f"{name}_test_auc"] = a
        res[f"{name}_test_p"] = p
        print(f"  {name}: AUC={a:.4f} p={p:.2g}", flush=True)

    for nm, S in (("armA", S_A), ("armB", S_B), ("armC", S_Cq)):
        ev(S, nm)

    # 配对 bootstrap ΔAUC(C−A)（逐对 win-rate，V8 机器同构）
    wa = np.array([np.mean(S_Cq[a, b] > S_Cq[a, negs[qgenes[a]]]) for a, b in pos_idx])
    wb = np.array([np.mean(S_A[a, b] > S_A[a, negs[qgenes[a]]]) for a, b in pos_idx])
    rb = np.random.default_rng(20260916)
    idx = np.arange(len(wa))
    ds = np.empty(4000)
    for t in range(4000):
        s = rb.choice(idx, len(idx), replace=True)
        ds[t] = wa[s].mean() - wb[s].mean()
    res["boot_C_vs_A"] = {"mean": float(ds.mean()),
                          "ci95": [float(np.percentile(ds, 2.5)),
                                   float(np.percentile(ds, 97.5))]}
    res["gate_primary"] = bool(res["armC_test_p"] < 0.05 and
                               res["boot_C_vs_A"]["ci95"][0] > 0)

    pd.DataFrame({"human": [h for h, _ in cand], "yeast": [y for _, y in cand],
                  "sources": ["|".join(sorted(pairs_src[(h, y)])) for h, y in cand],
                  "score_A": [S_A[qi[h], pool_pos[y]] for h, y in cand],
                  "score_C": [S_Cq[qi[h], pool_pos[y]] for h, y in cand]}
                 ).to_csv(OUT / "ev0_pairs.tsv", sep="\t", index=False)
    res["inputs_sha"] = {f: sha256_file(str(MAP / f"{f}.tsv"))
                         for f in ("orthodb_yeast_human", "oma_yeast_human",
                                   "inparanoid_yeast_human")}
    (OUT / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"GATE: {res['gate_primary']} | pairs={len(cand)} genes={len(qgenes)}",
          flush=True)


if __name__ == "__main__":
    main()
