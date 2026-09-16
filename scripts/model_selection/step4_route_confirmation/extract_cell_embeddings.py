#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_cell_embeddings.py — 路线A: 用微调后的 scGPT 提取 GSE125162 细胞 embedding(T1 输入)

用法:
    cd /public/home/mengxl/dzy/yeastbridge
    CUDA_VISIBLE_DEVICES=1 /public/home/mengxl/dzy/envs/yeastbridge/bin/python \
        scripts/routeA/extract_cell_embeddings.py            # 全量 38,225 细胞
    ... extract_cell_embeddings.py --max-cells 2000          # 冒烟

流程(全部官方组件):
    语料缓存 data/routeA/corpus_sc125162_norm.pt(normalize_total+log1p, 与微调一致)
    -> scgpt.preprocess.binning(n_bins=51) 分箱
    -> 每细胞按分箱值降序取 top-1200 基因(确定性截断, 对齐训练 max_seq_len; 并列由 topk 决定)
    -> 官方 TransformerModel.encode_batch(均值池化 + L2 归一化, 模型内建行为)
输出:
    data/routeA/cell_emb_gse125162.npz  {emb: (N,512) float32, obs_names: (N,)}
"""
import argparse
import inspect
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
ASSETS = Path("/public/home/mengxl/dzy/yeastbridge")  # legacy frozen asset root (scGPT weights/h5ad/mappings/external src); outputs go to this project
sys.path.insert(0, str(ASSETS / "src" / "external" / "scGPT"))

import anndata as ad  # noqa: E402
from scgpt.model import TransformerModel  # noqa: E402
from scgpt.preprocess import binning  # noqa: E402
from scgpt.tokenizer import GeneVocab  # noqa: E402

RA = ROOT / "results_model_selection" / "step4_route_confirmation" / "routeA_assets"
OUT = RA / "cell_emb_gse125162.npz"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--max-seq-len", type=int, default=1200)
    ap.add_argument("--max-cells", type=int, default=0)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--ckpt", default="routeA/scgpt_yeast_ft.pt",
                    help="models/ 下的权重相对路径(meta 同名 .meta.json);"
                         "如 routeA/scgpt_yeast_ft_holdout.pt、routeC/scgpt_yeast_scratch.pt")
    ap.add_argument("--out", default="results_model_selection/step4_route_confirmation/routeA_assets/cell_emb_gse125162.npz",
                    help="输出 npz 的项目相对路径")
    ap.add_argument("--protein-init", action="store_true",
                    help="路线B: 按蛋白注入层重建模型(scripts/routeB/protein_inject.py)")
    args = ap.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    vocab = GeneVocab.from_dict(json.load(open(RA / "routeA_vocab.json")), default_token="<pad>")
    ckpt_path = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / args.ckpt
    meta = json.load(open(ckpt_path.with_suffix(".meta.json")))["config_used"]
    sig = inspect.signature(TransformerModel.__init__).parameters
    kw = dict(ntoken=len(vocab), d_model=meta["embsize"], nhead=meta["nheads"], d_hid=meta["d_hid"],
              nlayers=meta["nlayers"], dropout=meta["dropout"], pad_token="<pad>",
              pad_value=meta["pad_value"], do_mvc=True, input_emb_style="continuous",
              n_input_bins=meta["n_bins"], use_fast_transformer=False, vocab=vocab)
    model = TransformerModel(**{k: v for k, v in kw.items() if k in sig})
    if args.protein_init:
        sys.path.insert(0, str(ASSETS / "scripts" / "routeB"))
        import protein_inject
        genes = [g for g, _ in sorted(((g, i) for g, i in vocab.get_stoi().items()
                                       if not g.startswith("<")), key=lambda x: x[1])]
        mat, n_missing = protein_inject.build_esm2_matrix(ROOT, genes)
        model.encoder.embedding = protein_inject.ProteinEmbeddingInjector(
            mat, len(vocab) - len(genes), meta["embsize"])
        print(f"[routeB] 注入层重建(ESM2 缺失 {n_missing}); 权重由 ckpt 全量加载")
    model.load_state_dict(torch.load(ckpt_path, map_location="cpu", weights_only=True))
    model.eval().to(device)
    print(f"[model] 微调权重加载完成, device={device}")

    d = torch.load(RA / "corpus_sc125162_norm.pt", map_location="cpu", weights_only=True)
    expr, gene_ids = d["expr"], d["gene_ids"]
    if args.max_cells:
        expr = expr[: args.max_cells]
    n = expr.shape[0]
    k = args.max_seq_len
    print(f"[data] {n} 细胞, 分箱 + top-{k} 截断 ...")
    t0 = time.time()

    src = torch.empty(n, k, dtype=torch.long)
    val = torch.empty(n, k, dtype=torch.float32)
    for i in range(n):
        b = binning(row=expr[i], n_bins=meta["n_bins"])
        v, ix = torch.topk(b, k)
        src[i] = gene_ids[ix]
        val[i] = v
        if i % 10000 == 0 and i:
            print(f"  {i}/{n} ({time.time()-t0:.0f}s)", flush=True)

    print(f"[encode] 前向取 cell_emb(官方定义: 均值池化+L2归一化), batch_size={args.batch_size} ...")
    embs = []
    with torch.no_grad():
        for i in range(0, n, args.batch_size):
            g = src[i:i + args.batch_size].to(device)
            v = val[i:i + args.batch_size].to(device)
            kpm = torch.zeros(g.shape[0], g.shape[1], dtype=torch.bool, device=device)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device.type == "cuda"):
                out = model(g, v, src_key_padding_mask=kpm, MVC=False)
            embs.append(out["cell_emb"].float().cpu())
            if (i // args.batch_size) % 200 == 0:
                print(f"  batch {i//args.batch_size}/{(n+args.batch_size-1)//args.batch_size} "
                      f"({time.time()-t0:.0f}s)", flush=True)
    emb = torch.cat(embs).numpy()
    obs = ad.read_h5ad(ASSETS / "data" / "processed" / "sc_gse125162" / "counts_raw.h5ad", backed="r").obs_names.to_numpy()
    if args.max_cells:
        obs = obs[: args.max_cells]
    assert len(obs) == len(emb)
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, emb=emb.astype(np.float32), obs_names=obs)
    print(f"[out] {out}  emb {emb.shape}  NaN: {bool(np.isnan(emb).any())}  ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
