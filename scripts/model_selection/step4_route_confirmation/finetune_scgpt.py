#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finetune_scgpt.py — 路线A: scGPT 换酵母词表 + 同源初始化 + 酵母语料继续预训练

用法:
    cd /public/home/mengxl/dzy/YeastBridge_0912
    CUDA_VISIBLE_DEVICES=1 /public/home/mengxl/dzy/YeastBridge_0912/env2/bin/python \
        scripts/routeA/finetune_scgpt.py            # 正式训练(后台)
    ... finetune_scgpt.py --smoke                 # 冒烟:512 细胞 60 步, 不记 mlflow

做法(官方源码优先, 全部组件来自 src/external/scGPT):
    1. 语料: GSE125162 counts_raw.h5ad (38,225 细胞 × 6,378 systematic 基因, 真实 counts)。
       预处理 = 官方 Preprocessor 流程: normalize_total(1e4) + log1p; 分箱交由官方
       DataCollator(do_binning, n_bins=51) 在 batch 内完成。缓存 data/routeA/corpus_sc125162_norm.pt。
       (SPELL 为 log2 ratio, 非 counts, 需专门的预处理, 本轮不混入——见 data/spell/README.md)
    2. 模型: 官方 TransformerModel, 结构参数全部取自 models/scgpt/scGPT_human/args.json
       (512d/12层/8头/d_hid 512/dropout 0.2/n_bins 51/pad_value -2/mask_value -1/MVC 开)。
       flash-attn 未装 -> use_fast_transformer=False(PyTorch 原生 attention, 结果等价)。
    3. 权重移植: best_model.pt 中所有形状匹配的键原样加载(Transformer 主干/value encoder/
       MVC decoder 等), encoder.embedding.weight 形状不同(60697->6736)跳过,
       改由 data/routeA/routeA_init_embeddings.npy 写入(2431 同源均值 + 4302 随机, seed=42)。
    4. 训练目标 = 官方生成式预训练: MLM masked-MSE + MVC masked-MSE(training_tasks=both);
       mask 比例每 batch 从 args.json 的 [0.25, 0.5, 0.75] 随机抽取;
       max_seq_len 1200, 超长随机采样(trunc_by_sample 官方行为)。
       与官方预训练的两处有意偏差: bf16 代替 fp16(A6000 上数值更稳); warmup 100 步代替 10000
       (本任务 ~7k 步, 原 warmup 比整个训练还长)。其余超参(lr 1e-4, batch 32, epochs 6)沿用 args.json。
    5. 产出:
       models/routeA/scgpt_yeast_ft.pt          最终权重(state_dict)
       models/routeA/scgpt_yeast_ft.meta.json   训练配置 + 权重移植对账
       results/routeA/finetune_report.json      loss 曲线摘要
       MLflow: experiment=yeastbridge_routeA, run=routeA_pretrain_scgpt_sc125162_s<seed>
"""
import argparse
import inspect
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[3]
ASSETS = Path("/public/home/mengxl/dzy/yeastbridge")  # legacy frozen asset root (scGPT weights/h5ad/mappings/external src); outputs go to this project
sys.path.insert(0, str(ASSETS / "src" / "external" / "scGPT"))

import anndata as ad  # noqa: E402
from scipy import sparse  # noqa: E402
from scgpt.data_collator import DataCollator  # noqa: E402
from scgpt.loss import masked_mse_loss  # noqa: E402
from scgpt.model import TransformerModel  # noqa: E402
from scgpt.tokenizer import GeneVocab  # noqa: E402
from scgpt.utils import set_seed  # noqa: E402

SCGPT_DIR = ASSETS / "models" / "scgpt" / "scGPT_human"
RA_DIR = ROOT / "results_model_selection" / "step4_route_confirmation" / "routeA_assets"
H5AD = ASSETS / "data" / "processed" / "sc_gse125162" / "counts_raw.h5ad"
CORPUS_CACHE = RA_DIR / "corpus_sc125162_norm.pt"
OUT_MODEL = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeA" / "scgpt_yeast_ft.pt"
OUT_META = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeA" / "scgpt_yeast_ft.meta.json"
OUT_REPORT = ROOT / "results_model_selection" / "step4_route_confirmation" / "ablations" / "finetune_report.json"
MLFLOW_URI = "sqlite:////public/home/mengxl/dzy/YeastBridge_0912/results_model_selection/step4_route_confirmation/mlflow.db"


def build_corpus(vocab, seed):
    """counts_raw.h5ad -> normalize_total(1e4)+log1p -> 稠密张量缓存。
    返回 (expr[N,G] float32, gene_ids[G] int64)。"""
    if CORPUS_CACHE.exists():
        print(f"[corpus] 加载缓存 {CORPUS_CACHE}")
        d = torch.load(CORPUS_CACHE, map_location="cpu", weights_only=True)
        return d["expr"], d["gene_ids"]
    print("[corpus] 首次构建: 读取 h5ad + 归一化 ...")
    t0 = time.time()
    a = ad.read_h5ad(H5AD)
    keep = [g in vocab for g in a.var.index]
    a = a[:, keep].copy()
    X = a.X if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    X = X.astype(np.float32).tocsr()
    rs = np.asarray(X.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    X = X.multiply(1e4 / rs[:, None]).tocsr()
    X.data = np.log1p(X.data)
    expr = torch.from_numpy(X.toarray().astype(np.float32))
    gene_ids = torch.tensor([vocab[g] for g in a.var.index], dtype=torch.long)
    CORPUS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"expr": expr, "gene_ids": gene_ids, "var_index": list(a.var.index)}, CORPUS_CACHE)
    print(f"[corpus] {expr.shape[0]} 细胞 × {expr.shape[1]} 基因, "
          f"{time.time()-t0:.0f}s, 缓存至 {CORPUS_CACHE}")
    return expr, gene_ids


class CorpusDataset(Dataset):
    """每例 {'id','genes','expressions'};expressions 克隆返回
    (官方 DataCollator 会原地写分箱结果, 不能污染缓存)。"""

    def __init__(self, expr, gene_ids):
        self.expr = expr
        self.gene_ids = gene_ids

    def __len__(self):
        return self.expr.shape[0]

    def __getitem__(self, idx):
        return {"id": idx, "genes": self.gene_ids, "expressions": self.expr[idx].clone()}


def build_model(vocab, cfg, device, from_scratch=False):
    sig = inspect.signature(TransformerModel.__init__).parameters
    kwargs = dict(
        ntoken=len(vocab), d_model=cfg["embsize"], nhead=cfg["nheads"],
        d_hid=cfg["d_hid"], nlayers=cfg["nlayers"], nlayers_cls=cfg["n_layers_cls"],
        dropout=cfg["dropout"], pad_token=cfg["pad_token"], pad_value=cfg["pad_value"],
        do_mvc=True, do_dab=False, input_emb_style=cfg["input_emb_style"],
        n_input_bins=cfg["n_bins"], use_fast_transformer=False, vocab=vocab,
    )
    kwargs = {k: v for k, v in kwargs.items() if k in sig}
    model = TransformerModel(**kwargs)

    if from_scratch:
        # 路线C: 构造器已执行官方 init_weights(embedding uniform(-0.1,0.1)), 不加载任何外部权重
        n_params = sum(p.numel() for p in model.parameters())
        report = {"from_scratch": True, "transferred_keys": 0, "init": "official init_weights"}
        print(f"[weights] 路线C 从零初始化(官方 init_weights), 0 键移植")
        return model.to(device), report

    sd_human = torch.load(SCGPT_DIR / "best_model.pt", map_location="cpu", weights_only=True)
    sd_self = model.state_dict()
    transferred = {k: v for k, v in sd_human.items()
                   if k in sd_self and sd_self[k].shape == v.shape
                   and k != "encoder.embedding.weight"}
    # 关键转换: 人类 ckpt 是 flash-attn 命名(Wqkv.*), 原生模型是 in_proj_*;
    # 两者同为 [q;k;v] 行堆叠 head-major 布局, 直接改名拷贝即可,
    # 数值等价性见 scripts/routeA/verify_wqkv_conversion.py (PASS)。
    n_wqkv = 0
    for i in range(cfg["nlayers"]):
        for src_suf, dst_suf in (("Wqkv.weight", "in_proj_weight"), ("Wqkv.bias", "in_proj_bias")):
            sk = f"transformer_encoder.layers.{i}.self_attn.{src_suf}"
            dk = f"transformer_encoder.layers.{i}.self_attn.{dst_suf}"
            if sk in sd_human and dk in sd_self and dk not in transferred:
                transferred[dk] = sd_human[sk]
                n_wqkv += 1
    missing = [k for k in sd_self if k not in transferred]
    model.load_state_dict(transferred, strict=False)
    init_emb = torch.from_numpy(np.load(RA_DIR / "routeA_init_embeddings.npy"))
    assert init_emb.shape == sd_self["encoder.embedding.weight"].shape
    model.encoder.embedding.weight.data.copy_(init_emb)

    n_core = sum(1 for k in transferred if k.startswith("transformer_encoder"))
    report = {
        "transferred_keys": len(transferred), "core_transformer_keys": n_core,
        "wqkv_to_inproj_converted": n_wqkv,
        "model_total_keys": len(sd_self),
        "human_ckpt_keys_dropped": [k for k in sd_human if k not in transferred],
        "fresh_keys_head": missing[:10], "n_fresh_keys": len(missing),
    }
    print(f"[weights] 移植 {len(transferred)}/{len(sd_self)} 键(主干 {n_core} 个, "
          f"含 Wqkv->in_proj 转换 {n_wqkv} 个); "
          f"新初始化 {len(missing)} 键; 人类 ckpt 弃用 {len(report['human_ckpt_keys_dropped'])} 键")
    return model.to(device), report


def run_epoch(model, loader, collators, vocab, cfg, device, rng,
              optimizer=None, scheduler=None, max_steps=None, tag="train"):
    training = optimizer is not None
    model.train() if training else model.eval()
    pad_id = vocab[cfg["pad_token"]]
    tot_loss = tot_mlm = tot_mvc = nb = 0
    for step, examples in enumerate(loader):
        if max_steps and step >= max_steps:
            break
        collator = collators[rng.integers(len(collators))] if training else collators[0]
        data = collator(examples)
        gene = data["gene"].to(device)
        target = data["expr"].to(device)
        masked = data["masked_expr"].to(device)
        mask = masked.eq(cfg["mask_value"])
        kpm = gene.eq(pad_id)
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device.type == "cuda"):
            out = model(gene, masked, src_key_padding_mask=kpm, MVC=True)
            loss_mlm = masked_mse_loss(out["mlm_output"].float(), target, mask)
            loss_mvc = masked_mse_loss(out["mvc_output"].float(), target, mask)
            loss = loss_mlm + loss_mvc
        if training:
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            if scheduler is not None:
                scheduler.step()
        tot_loss += loss.item(); tot_mlm += loss_mlm.item(); tot_mvc += loss_mvc.item(); nb += 1
        if training and step % 100 == 0:
            print(f"  [{tag}] step {step}  loss {loss.item():.4f} "
                  f"(mlm {loss_mlm.item():.4f} mvc {loss_mvc.item():.4f})", flush=True)
    return {"loss": tot_loss / max(nb, 1), "mlm": tot_mlm / max(nb, 1), "mvc": tot_mvc / max(nb, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--epochs", type=int, default=6)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--max-seq-len", type=int, default=1200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--warmup", type=int, default=100)
    ap.add_argument("--n-valid", type=int, default=500)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--holdout-conditions", default="",
                    help="逗号分隔的条件名, 这些条件的细胞完全不进训练/验证(T1 零样本公平性臂);"
                         "模型存为 scgpt_yeast_ft_holdout.pt 并写 data/routeA/holdout_conditions.json")
    ap.add_argument("--from-scratch", action="store_true",
                    help="路线C: 不加载任何人类权重, 官方 init_weights 随机初始化从零训练;"
                         "产出归 models/routeC/ 与 results/routeC/")
    ap.add_argument("--protein-init", action="store_true",
                    help="路线B: 主干移植人类 scGPT(同路线A), 基因 embedding 层换为 ESM2 蛋白注入层"
                         "(scripts/routeB/protein_inject.py);产出归 models/routeB/ 与 results/routeB/")
    args = ap.parse_args()
    assert not (args.from_scratch and (args.holdout_conditions or args.protein_init)), \
        "from_scratch 为独立路线, 不与 holdout/protein_init 同用"

    cfg = json.load(open(SCGPT_DIR / "args.json"))
    set_seed(args.seed)
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    vocab = GeneVocab.from_dict(json.load(open(RA_DIR / "routeA_vocab.json")),
                                default_token=cfg["pad_token"])
    expr, gene_ids = build_corpus(vocab, args.seed)
    hold = set()
    if args.holdout_conditions:
        hold = set(args.holdout_conditions.split(","))
        conds = ad.read_h5ad(H5AD, backed="r").obs["condition"].astype(str).to_numpy()
        assert len(conds) == len(expr), "obs 与语料缓存行数不一致"
        keep = np.array([c not in hold for c in conds])
        print(f"[holdout] 留出 {sorted(hold)}: 排除 {(~keep).sum()} 细胞, 余 {keep.sum()}")
        expr = expr[torch.from_numpy(np.nonzero(keep)[0])]
        RA_DIR.joinpath("holdout_conditions.json").write_text(json.dumps({
            "conditions": sorted(hold), "n_cells_excluded": int((~keep).sum())}, indent=1))
    if args.smoke:
        expr = expr[:512]

    # 固定划分: 种子置换后取末尾 n_valid 个细胞做验证(h5ad 行序按条件聚簇, 不能顺序切)
    n_valid = 0 if args.smoke else min(args.n_valid, len(expr) // 10)
    perm = np.random.default_rng(args.seed).permutation(len(expr))
    tr_idx, va_idx = perm[: len(expr) - n_valid], perm[len(expr) - n_valid:]
    train_ds = CorpusDataset(expr[tr_idx], gene_ids)
    valid_ds = CorpusDataset(expr[va_idx], gene_ids) if n_valid else None
    print(f"[data] train {len(train_ds)}  valid {n_valid}  genes {len(gene_ids)}")

    mk = lambda p: DataCollator(  # noqa: E731
        do_padding=True, pad_token_id=vocab[cfg["pad_token"]], pad_value=cfg["pad_value"],
        do_mlm=True, do_binning=True, mlm_probability=p, mask_value=cfg["mask_value"],
        max_length=args.max_seq_len, sampling=True, keep_first_n_tokens=0)
    collators = [mk(p) for p in cfg["mask_ratio"]]          # 训练: 每 batch 随机抽比例
    valid_collators = [mk(0.25)]                            # 验证: 固定 0.25 保证可比
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              collate_fn=lambda ex: ex, num_workers=0, drop_last=True)
    valid_loader = DataLoader(valid_ds, batch_size=64, shuffle=False,
                              collate_fn=lambda ex: ex, num_workers=0) if valid_ds else None

    if args.protein_init:
        sys.path.insert(0, str(ASSETS / "scripts" / "routeB"))
        import protein_inject
        model, wrep = protein_inject.build_protein_model(vocab, cfg, device, ROOT, build_model)
    else:
        model, wrep = build_model(vocab, cfg, device, from_scratch=args.from_scratch)
    n_params = sum(p.numel() for p in model.parameters())
    print(f"[model] {n_params/1e6:.1f}M 参数, device={device}")

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lambda s: min(1.0, (s + 1) / args.warmup))

    epochs = 1 if args.smoke else args.epochs
    max_steps = 60 if args.smoke else None
    rng = np.random.default_rng(args.seed)
    history = []
    t0 = time.time()
    for ep in range(epochs):
        tr = run_epoch(model, train_loader, collators, vocab, cfg, device, rng,
                       optimizer, sched, max_steps, tag=f"ep{ep}")
        row = {"epoch": ep, **{f"train_{k}": v for k, v in tr.items()}}
        if valid_loader:
            with torch.no_grad():
                va = run_epoch(model, valid_loader, valid_collators, vocab, cfg, device, rng)
            row.update({f"valid_{k}": v for k, v in va.items()})
        history.append(row)
        print(f"[epoch {ep}] " + "  ".join(f"{k}={v:.4f}" for k, v in row.items() if k != "epoch")
              + f"  ({time.time()-t0:.0f}s)", flush=True)

    report = {"args": vars(args), "weight_transfer": wrep, "n_params": n_params,
              "history": history, "wall_time_s": round(time.time() - t0, 1)}
    if args.smoke:
        print("[smoke] 完成, 不落盘不记 mlflow")
        print(json.dumps(history, indent=1))
        return 0

    if args.from_scratch:
        suf = ""
        out_model = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeC" / "scgpt_yeast_scratch.pt"
        out_meta = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeC" / "scgpt_yeast_scratch.meta.json"
        out_report = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeC" / "pretrain_report.json"
        experiment = "yeastbridge_routeC"
    elif args.protein_init:
        suf = "_holdout" if hold else ""
        out_model = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeB" / f"scgpt_yeast_protein{suf}.pt"
        out_meta = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeB" / f"scgpt_yeast_protein{suf}.meta.json"
        out_report = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeB" / f"finetune_report{suf}.json"
        experiment = "yeastbridge_routeB"
    else:
        suf = "_holdout" if hold else ""
        out_model = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeA" / f"scgpt_yeast_ft{suf}.pt"
        out_meta = ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeA" / f"scgpt_yeast_ft{suf}.meta.json"
        out_report = ROOT / "results_model_selection" / "step4_route_confirmation" / "ablations" / f"finetune_report{suf}.json"
        experiment = "yeastbridge_routeA"
    out_model.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), out_model)
    out_meta.write_text(json.dumps({"config_used": {
        "ntoken": len(vocab), "embsize": cfg["embsize"], "nheads": cfg["nheads"],
        "d_hid": cfg["d_hid"], "nlayers": cfg["nlayers"], "dropout": cfg["dropout"],
        "n_bins": cfg["n_bins"], "pad_token": cfg["pad_token"], "pad_value": cfg["pad_value"],
        "mask_value": cfg["mask_value"], "input_emb_style": cfg["input_emb_style"],
        "do_mvc": True, "use_fast_transformer": False, "max_seq_len": args.max_seq_len,
    }, "weight_transfer": wrep, "train_args": vars(args),
        "holdout_conditions": sorted(hold)}, indent=1))
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(json.dumps(report, indent=1))

    if args.protein_init:
        # 训练后基因 embedding 表(注入层展开), 供 eval 的 routeb_protein 特征直接使用
        gt = protein_inject.gene_embedding_table(model)
        np.save(ROOT / "results_model_selection" / "step4_route_confirmation" / "models" / "routeB" / "routeb_gene_embeddings.npy", gt)
        print(f"[out] {ROOT}/models/routeB/routeb_gene_embeddings.npy {gt.shape}")

    import mlflow
    mlflow.set_tracking_uri(MLFLOW_URI)
    mlflow.set_experiment(experiment)
    route_tag = "routeC" if args.from_scratch else ("routeB" if args.protein_init else "routeA")
    run_name = f"{route_tag}_pretrain_scgpt_sc125162_s{args.seed}{suf}"
    with mlflow.start_run(run_name=run_name):
        mlflow.log_params({**{k: v for k, v in vars(args).items() if k != "smoke"},
                           "n_params": n_params, "corpus": "gse125162_counts",
                           "transferred_keys": wrep["transferred_keys"]})
        for row in history:
            mlflow.log_metrics({k: v for k, v in row.items() if k != "epoch"}, step=row["epoch"])
        mlflow.log_artifact(str(out_meta))
        mlflow.log_artifact(str(out_report))
    print(f"[out] {out_model}\n[out] {out_meta}\n[out] {out_report}\n[mlflow] {run_name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
