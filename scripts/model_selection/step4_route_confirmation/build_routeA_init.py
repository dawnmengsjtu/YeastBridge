#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_routeA_init.py — 路线A 酵母基因 embedding 初始化表构建

用法:
    cd /public/home/mengxl/dzy/YeastBridge_0912
    /public/home/mengxl/dzy/YeastBridge_0912/env2/bin/python scripts/routeA/build_routeA_init.py

输入:
    data/mappings/gene_master.tsv                  酵母基因总表(systematic 为主键, 决定输出行序)
    data/mappings/orthodb_yeast_human_s288c.tsv    (yeast_sgd, human_symbol)
    data/mappings/oma_yeast_human.tsv              (yeast_systematic, human_symbol, ortholog_type)
    data/mappings/inparanoid_yeast_human.tsv       (yeast_systematic, human_symbol, ortholog_type)
    models/scgpt/scGPT_human/vocab.json            scGPT 人类全基因组 vocab(符号 -> 行号)
    models/scgpt/scGPT_human/best_model.pt         scGPT 预训练权重(取 encoder.embedding.weight)

初始化规则:
    1. 有同源且至少一个人类符号在 scGPT vocab 内 -> 取全部命中符号的 embedding 均值(init_type=ortholog)
    2. 其余(无同源/符号不在 vocab)-> 按人类 embedding 矩阵逐维 mean/std 高斯采样(init_type=random, seed=42)
    3. 特殊符号 <pad>/<cls>/<eoc>(若人类 vocab 中存在)直接拷贝对应 embedding, 追加在基因之后

输出(data/routeA/):
    routeA_vocab.json           酵母基因 vocab(systematic -> 行号), 特殊符号追加在末尾
    routeA_init_embeddings.npy  float32 (len(vocab), 512), 行序与 routeA_vocab.json 对齐
    routeA_init.tsv             每基因元数据(init_type, 命中符号, 来源, 是否1:1)
    README.md                   构建方法 + 统计 + 复现命令
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
ASSETS = Path("/public/home/mengxl/dzy/yeastbridge")  # legacy frozen asset root (scGPT weights/h5ad/mappings/external src); outputs go to this project
M = ASSETS / "data" / "mappings"
SCGPT = ASSETS / "models" / "scgpt" / "scGPT_human"
OUT = ROOT / "results_model_selection" / "step4_route_confirmation" / "routeA_assets"

SEED = 42
SPECIALS = ["<pad>", "<cls>", "<eoc>"]

SOURCES = [
    ("OrthoDB", M / "orthodb_yeast_human_s288c.tsv", "yeast_sgd", "human_symbol", None),
    ("OMA", M / "oma_yeast_human.tsv", "yeast_systematic", "human_symbol", "ortholog_type"),
    ("InParanoid", M / "inparanoid_yeast_human.tsv", "yeast_systematic", "human_symbol", "ortholog_type"),
]


def load_master():
    genes = []
    with open(M / "gene_master.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = row["systematic"].strip()
            if g:
                genes.append((g, row["common"].strip()))
    return genes


def load_orthologs():
    """yeast systematic -> {symbols: set, sources: set, one_to_one: bool}"""
    table = {}
    for name, path, yc, hc, tc in SOURCES:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                y, h = row[yc].strip(), row[hc].strip()
                if not y or not h:
                    continue
                rec = table.setdefault(y, {"symbols": set(), "sources": set(), "one_to_one": False})
                rec["symbols"].add(h)
                rec["sources"].add(name)
                if tc and row.get(tc, "").strip() == "1:1":
                    rec["one_to_one"] = True
    return table


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    master = load_master()
    ortho = load_orthologs()

    vocab_human = json.load(open(SCGPT / "vocab.json"))
    sd = torch.load(SCGPT / "best_model.pt", map_location="cpu", weights_only=True)
    emb = sd["encoder.embedding.weight"].numpy().astype(np.float32)
    assert emb.shape[0] == len(vocab_human), "vocab 与 embedding 行数不一致"
    emb_dim = emb.shape[1]
    print(f"[load] master genes={len(master)}  ortholog yeast genes={len(ortho)}")
    print(f"[load] scGPT human vocab={len(vocab_human)}  emb_dim={emb_dim}")

    dim_mean = emb.mean(axis=0)
    dim_std = emb.std(axis=0)
    rng = np.random.default_rng(SEED)

    n_genes = len(master)
    out_emb = np.empty((n_genes, emb_dim), dtype=np.float32)
    meta_rows = []
    n_ortholog = n_random = 0

    for i, (sysname, common) in enumerate(master):
        rec = ortho.get(sysname)
        used = sorted(s for s in rec["symbols"] if s in vocab_human) if rec else []
        if used:
            out_emb[i] = emb[[vocab_human[s] for s in used]].mean(axis=0)
            init_type = "ortholog"
            n_ortholog += 1
        else:
            out_emb[i] = rng.normal(dim_mean, dim_std).astype(np.float32)
            init_type = "random"
            n_random += 1
        meta_rows.append({
            "systematic": sysname,
            "common": common,
            "init_type": init_type,
            "n_ortholog_symbols": len(rec["symbols"]) if rec else 0,
            "n_symbols_in_vocab": len(used),
            "human_symbols_used": ",".join(used),
            "sources": ",".join(sorted(rec["sources"])) if rec else "",
            "one_to_one": int(rec["one_to_one"]) if rec else 0,
        })

    # 特殊符号: 直接拷贝人类 embedding, 追加在基因之后
    vocab_out = {g: i for i, (g, _) in enumerate(master)}
    extra = []
    for sp in SPECIALS:
        if sp in vocab_human:
            vocab_out[sp] = len(vocab_out)
            extra.append(emb[vocab_human[sp]])
    if extra:
        out_emb = np.vstack([out_emb, np.stack(extra).astype(np.float32)])

    np.save(OUT / "routeA_init_embeddings.npy", out_emb)
    with open(OUT / "routeA_vocab.json", "w") as f:
        json.dump(vocab_out, f, indent=0)
    with open(OUT / "routeA_init.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(meta_rows[0].keys()))
        w.writeheader()
        w.writerows(meta_rows)

    # ---- 报告 ----
    n_one2one = sum(r["one_to_one"] and r["init_type"] == "ortholog" for r in meta_rows)
    pre_norm = np.linalg.norm(emb, axis=1).mean()
    ortho_idx = [i for i, r in enumerate(meta_rows) if r["init_type"] == "ortholog"]
    rand_idx = [i for i, r in enumerate(meta_rows) if r["init_type"] == "random"]
    print(f"\n[result] ortholog-init: {n_ortholog} ({n_ortholog/n_genes:.1%})  其中含1:1同源: {n_one2one}")
    print(f"[result] random-init:   {n_random} ({n_random/n_genes:.1%})")
    print(f"[check] 平均向量模长: 预训练={pre_norm:.3f}  "
          f"ortholog={np.linalg.norm(out_emb[ortho_idx], axis=1).mean():.3f}  "
          f"random={np.linalg.norm(out_emb[rand_idx], axis=1).mean():.3f}")

    # 点检: CDC28 (YBR160W) 应当映射到 CDK1 家族
    for r in meta_rows:
        if r["systematic"] == "YBR160W":
            print(f"[spot] YBR160W/CDC28 -> {r['human_symbols_used']} ({r['init_type']}, sources={r['sources']})")
            i = vocab_out["YBR160W"]
            if r["human_symbols_used"] == "CDK1":
                exact = np.allclose(out_emb[i], emb[vocab_human["CDK1"]], atol=1e-6)
                print(f"[spot] 单一同源时与人类 CDK1 行完全一致: {exact}")
            break

    readme = f"""# 路线A 酵母基因 embedding 初始化表

构建: `python scripts/routeA/build_routeA_init.py` (seed={SEED})

- 基因全集: data/mappings/gene_master.tsv ({n_genes} 个, 行序即 vocab 行序)
- 同源并集: OrthoDB + OMA + InParanoid, 人类符号命中 scGPT vocab 者取均值初始化
- ortholog-init: {n_ortholog} ({n_ortholog/n_genes:.1%});  random-init: {n_random} ({n_random/n_genes:.1%})
- random-init: 按人类 embedding 逐维 mean/std 高斯采样(np.random.default_rng({SEED}))
- 特殊符号 {', '.join(sp for sp in SPECIALS if sp in vocab_out)} 从人类 vocab 直接拷贝, 追加在基因之后
- 文件: routeA_vocab.json / routeA_init_embeddings.npy (float32, {out_emb.shape}) / routeA_init.tsv
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(f"\n[out] {OUT}/routeA_vocab.json, routeA_init_embeddings.npy {out_emb.shape}, routeA_init.tsv, README.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
