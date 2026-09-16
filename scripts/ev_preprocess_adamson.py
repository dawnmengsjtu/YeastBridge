#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EV1 预处理：Adamson GSE90546 (3 批次恢复件) -> 外样本查询集。

口径（EV_PROTOCOL.md §EV1，先冻结后运行）：
  - 输入：GSE90546_RAW.tar（3 批次 10x：10X001/10X005/10X010）+ 同目录 sidecar identities。
  - 每批次：good coverage TRUE 细胞按 guide identity 聚 pseudobulk；
    guide 目标 = identity 下划线前段；非基因符号/对照（含 '62(mod)' 等）入对照池。
  - 签名 = 各批次 log1p(CP10K) 目标池 - 对照池，逐基因；跨批次取均值（记录批次数）。
  - 标签 = Norman 207 单基因签名的 gene->function_label 查表直传（冻结表，运行时 sha 记录）；
    查表未命中的基因保留在矩阵中但 metadata 标 label_unmapped=true，不进 EV1 门端点。
  - split_group = 基因符号（LOGO 单位）；split 留空（外部集整集为 test 域）。
  - 输出 sha256 全记录；矩阵/元数据列名对齐 cross_species_match.load_matrix。
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

EV = Path("/public/home/mengxl/dzy/YeastBridge_0912")
RAW = Path("/public/home/mengxl/dzy/YeastBridge_0912/raw/"
           "externalvalidation/adamson_raw")
EXT = EV / "raw/externalvalidation/adamson_extracted"
OUT = EV / "results_a6000/externalvalidation/ev1_adamson"
OUT.mkdir(parents=True, exist_ok=True)
NORMAN_META = Path("/public/home/mengxl/dzy/YeastBridge_0912/results_a6000")  # 占位，下方解析
BATCHES = ["GSM2406675_10X001", "GSM2406677_10X005", "GSM2406681_10X010"]
MIN_CELLS_PER_TARGET = 5


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def load_batch(tag):
    genes = pd.read_csv(EXT / f"{tag}_genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    barcodes = [l.strip() for l in gzip.open(EXT / f"{tag}_barcodes.tsv.gz", "rt")]
    m = mmread(str(EXT / f"{tag}_matrix.mtx.txt.gz")).tocsr()
    # 10x mtx 为 genes x cells；统一转为 cells x genes
    if m.shape[0] == len(genes) and m.shape[1] == len(barcodes):
        m = m.T.tocsr()
    assert m.shape[0] == len(barcodes) and m.shape[1] == len(genes), \
        f"mtx 形状不符: {tag} {m.shape} bc={len(barcodes)} g={len(genes)}"
    ident = pd.read_csv(RAW / f"{tag}_cell_identities.csv.gz")
    ident = ident[ident["good coverage"].astype(str).str.upper() == "TRUE"]
    bc2guide = dict(zip(ident["cell BC"], ident["guide identity"].astype(str)))
    return genes, barcodes, m, bc2guide


def main():
    report = {"mode": "ev1_preprocess_adamson", "batches": BATCHES,
              "min_cells_per_target": MIN_CELLS_PER_TARGET,
              "inputs": {}}
    # tar + identities + genes 哈希
    report["inputs"]["tar"] = {"path": str(RAW / "GSE90546_RAW.tar"),
                               "sha256": sha(RAW / "GSE90546_RAW.tar")}
    for tag in BATCHES:
        report["inputs"][f"{tag}_identities"] = {
            "path": str(RAW / f"{tag}_cell_identities.csv.gz"),
            "sha256": sha(RAW / f"{tag}_cell_identities.csv.gz")}

    # ---- Norman 标签冻结表（查表直传；路径与 sha 对冻结配置断言）----
    import json as _json
    cfg51 = _json.loads(Path(
        "/public/home/mengxl/dzy/YeastBridge_0912/configs_a6000_frozen/"
        "v5_1_joint_formal.json").read_text())
    stage_root = Path("/public/home/mengxl/dzy/YeastBridge_0912/raw/externalvalidation/frozen_inputs")
    rel = cfg51["inputs"]["higher_metadata"]["path"]
    label_src = stage_root / rel
    if not label_src.is_file():  # A6000 兜底：stage 根即 mvp stage 目录
        label_src = Path(
            "/public/home/mengxl/dzy/YeastBridge_0912/raw/externalvalidation/frozen_inputs") / rel
    got_sha = sha(label_src)
    want_sha = cfg51["inputs"]["higher_metadata"]["sha256"]
    assert got_sha == want_sha, (
        f"Norman metadata sha 不符冻结配置: got {got_sha[:12]} want {want_sha[:12]}")
    nm = pd.read_csv(label_src, sep="\t")
    nm = nm[nm.benchmark_eligible.astype(str).str.lower() == "true"]
    single = nm[nm.perturbation_genes.astype(str).str.match(r"^[A-Za-z0-9@.-]+$")]
    gene2label = dict(zip(single.perturbation_genes.astype(str),
                          single.function_label.astype(str)))
    report["label_table"] = {"source": str(label_src),
                             "sha256": sha(label_src),
                             "n_signatures": int(len(nm)),
                             "n_single_gene": int(len(single)),
                             "n_unique_genes": len(gene2label)}
    print(f"label table frozen: {len(gene2label)} genes from {len(single)} sigs",
          flush=True)

    # ---- 逐批次 pseudobulk（基因空间按并集对齐，缺基因记 NaN）----
    deltas = {}       # target -> list of union-space delta (含 NaN)
    gene_union, gene_seen = [], set()
    batch_gene_space = {}
    per_batch_stats = {}
    for tag in BATCHES:
        genes, barcodes, m, bc2guide = load_batch(tag)
        syms = genes.symbol.tolist()
        batch_gene_space[tag] = syms
        for g in syms:
            if g and g not in gene_seen:
                gene_seen.add(g)
                gene_union.append(g)
        print(f"{tag}: cells={m.shape[0]}, genes={len(syms)}, "
              f"union={len(gene_union)}", flush=True)
        bc_index = {b: i for i, b in enumerate(barcodes)}
        rows, target_cells = [], {}
        ctrl_rows = []
        for bc, guide in bc2guide.items():
            i = bc_index.get(bc)
            if i is None:
                continue
            tgt = guide.split("_")[0]
            # 对照判定：非人类基因符号形态（数字开头/含括号/小写起始等）
            if not tgt[:1].isalpha() or tgt[:1].islower() or "(" in tgt:
                ctrl_rows.append(i)
            else:
                target_cells.setdefault(tgt, []).append(i)
        def pb_sum(idxs):
            if not idxs:
                return None
            return np.asarray(m[idxs].sum(axis=0)).ravel()
        ctrl = pb_sum(ctrl_rows)
        stats = {"cells_total": m.shape[0], "cells_used": len(rows) or
                 sum(len(v) for v in target_cells.values()) + len(ctrl_rows),
                 "control_cells": len(ctrl_rows), "targets": {}}
        n_target_kept = 0
        u_map = {g: i for i, g in enumerate(gene_union)}
        u_idx = np.array([u_map.get(g, -1) for g in syms])
        for tgt, idxs in sorted(target_cells.items()):
            if len(idxs) < MIN_CELLS_PER_TARGET:
                stats["targets"][tgt] = {"cells": len(idxs), "kept": False}
                continue
            s = pb_sum(idxs)
            if s is None or ctrl is None:
                continue
            lib_t, lib_c = s.sum(), ctrl.sum()
            if lib_t <= 0 or lib_c <= 0:
                continue
            d_batch = (np.log1p(s / lib_t * 1e4) - np.log1p(ctrl / lib_c * 1e4))
            d = np.full(len(gene_union), np.nan)
            ok = u_idx >= 0
            d[u_idx[ok]] = d_batch[ok]
            deltas.setdefault(tgt, []).append(d)
            stats["targets"][tgt] = {"cells": len(idxs), "kept": True}
            n_target_kept += 1
        stats["n_targets_kept"] = n_target_kept
        per_batch_stats[tag] = stats
        print(f"{tag}: cells={m.shape[0]}, ctrl={len(ctrl_rows)}, "
              f"targets_kept={n_target_kept}", flush=True)

    # ---- 跨批次 nan 均值 + 输出 ----
    targets = sorted(deltas)
    G = len(gene_union)
    X = np.zeros((len(targets), G))
    nb = np.zeros(len(targets), dtype=int)
    for i, t in enumerate(targets):
        stack = np.vstack(deltas[t])
        X[i] = np.nan_to_num(np.nanmean(stack, axis=0))
        nb[i] = stack.shape[0]
    # 列去重（同符号基因取首个）
    keep_col, seen = [], set()
    for j, g in enumerate(gene_union):
        if g and g not in seen:
            seen.add(g)
            keep_col.append(j)
    X = X[:, keep_col]
    cols = [gene_union[j] for j in keep_col]
    meta = pd.DataFrame({
        "signature_id": [f"adamson_{t}" for t in targets],
        "species": "human",
        "perturbation_genes": targets,
        "n_batches": nb,
        "label_unmapped": [t not in gene2label for t in targets],
        "function_label": [gene2label.get(t, "") for t in targets],
        "split_group": [f"adamson_{t}" for t in targets],
        "split": "",
        "assay_type": "Perturb-seq CRISPRi",
        "source": "GSE90546",
    })
    Xz = ((X - X.mean(axis=1, keepdims=True)) /
          np.where(X.std(axis=1, keepdims=True) > 0, X.std(axis=1, keepdims=True), 1))
    mat = pd.DataFrame(Xz, columns=cols)
    mat.insert(0, "signature_id", meta.signature_id.values)
    mat.to_csv(OUT / "adamson_matrix.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "adamson_metadata.tsv", sep="\t", index=False)
    report["per_batch"] = per_batch_stats
    report["outputs"] = {
        "matrix": {"path": str(OUT / "adamson_matrix.tsv"), "shape": list(Xz.shape)},
        "metadata": {"path": str(OUT / "adamson_metadata.tsv"),
                     "n": len(meta),
                     "n_label_mapped": int((~meta.label_unmapped).sum())},
    }
    (OUT / "preprocess_report.json").write_text(json.dumps(report, indent=1))
    print(f"[done] {len(meta)} targets, mapped {int((~meta.label_unmapped).sum())} "
          f"-> {OUT}", flush=True)


if __name__ == "__main__":
    main()
