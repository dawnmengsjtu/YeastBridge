#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Portable B2 training entry adapted from the recovered original script.

The original recipe remains available through implementation=historical-v1.
The submission-v2 recipe excludes per-row padding/special tokens from masking
and evaluates fixed validation masks in eval mode without gradients.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[1]
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("--config", type=Path, default=ROOT / "configs/b2_training.json")
_pre_args, _ = _pre.parse_known_args()
CFG = json.loads(_pre_args.config.read_text())
YB = (ROOT / CFG["yeastbridge_root"]).resolve()
ASSETS = ROOT / CFG["out_assets_dir"]
OUTDIR = ROOT / CFG["out_routes_dir"]
MODEL_DIR = YB / CFG["scf_model_dir"]
sys.path.insert(0, str(MODEL_DIR))

from load import gatherData, load_model_frommmf  # noqa: E402

N_GENES = CFG[
    "n_genes"
]  # gene_master 真实基因数(旧 routeA_vocab 6736 = 6733 基因 + 3 特殊符号)
SEQ_LEN = N_GENES + 2  # + 2 分辨率位; pos_emb 行数 = SEQ_LEN + 1(pad 位)
MAX_SEQ = CFG[
    "max_seq"
]  # 原 finetune_scgpt.py train_args.max_seq_len=1200, trunc_by_sample 官方行为
EPOCHS = CFG["epochs"]
BATCH = CFG["batch"]
ACCUM = CFG["accum"]
LR = CFG["lr"]
CLIP = CFG["clip"]
SEED = CFG["seed"]
MASK_P = CFG["mask_p"]
ZERO_MASK_P = CFG["zero_mask_p"]
FIXED = CFG.get("implementation", "historical-v1") == "submission-v2"


class ProteinEmbeddingInjector(nn.Module):
    """逐字复制自 scripts/routeB/protein_inject.py,仅 d_model 由调用方给定。"""

    def __init__(self, prot_mat: np.ndarray, n_special: int, d_model: int):
        super().__init__()
        self.n_genes = prot_mat.shape[0]
        self.register_buffer("prot", torch.from_numpy(prot_mat.astype(np.float32)))
        self.proj = nn.Linear(prot_mat.shape[1], d_model)
        self.special = nn.Embedding(n_special, d_model)

    def forward(self, idx: torch.Tensor) -> torch.Tensor:
        is_special = idx >= self.n_genes
        gene_idx = idx.clamp(max=self.n_genes - 1)
        out = self.proj(self.prot[gene_idx])
        if bool(is_special.any()):
            sp = self.special((idx - self.n_genes).clamp(min=0))
            out = torch.where(is_special.unsqueeze(-1), sp, out)
        return out


class Checkpointed(nn.Module):
    """梯度检查点包装:训练时重算换显存(否则 token_emb 可训练会把 12 层注意力
    矩阵全部驻留,~47GB)。use_reentrant=False 下末 2 层可训练参数照常收梯度。"""

    def __init__(self, m):
        super().__init__()
        self.m = m

    def forward(self, x, padding_mask=None):
        if self.training and x.requires_grad:
            return torch.utils.checkpoint.checkpoint(
                self.m, x, padding_mask=padding_mask, use_reentrant=False
            )
        return self.m(x, padding_mask=padding_mask)


def encdec_subset(data, data_raw, gene_ids, config):
    """load.getEncoerDecoderData 的逐字逻辑,唯一改动:data_gene_ids 由调用方给出
    (子集行的真实基因列号),不再假设 arange(全词表宽)。"""
    decoder_data = data.clone().detach()
    decoder_data_padding = (
        (gene_ids == config["seq_len"])
        if FIXED
        else torch.full_like(data, False, dtype=torch.bool)
    )
    encoder_data_labels = data_raw > 0
    encoder_data, encoder_data_padding = gatherData(
        decoder_data, encoder_data_labels, config["pad_token_id"]
    )
    new_data_raw = data_raw
    encoder_position_gene_ids, _ = gatherData(
        gene_ids, encoder_data_labels, config["pad_token_id"]
    )
    decoder_position_gene_ids = gene_ids.clone()
    data_mask_labels = None
    encoder_position_gene_ids[encoder_data_padding] = config["seq_len"]
    decoder_position_gene_ids[decoder_data_padding] = config["seq_len"]
    return (
        encoder_data,
        encoder_position_gene_ids,
        encoder_data_padding,
        encoder_data_labels,
        decoder_data,
        decoder_data_padding,
        new_data_raw,
        data_mask_labels,
        decoder_position_gene_ids,
    )


def pick_device(requested):
    if requested == "auto":
        if torch.cuda.is_available() and torch.cuda.mem_get_info()[0] / 1024**3 >= 8.0:
            return torch.device("cuda")
        return torch.device("cpu")
    return torch.device(requested)


def build_model(route, device):
    model, cfg = load_model_frommmf(str(YB / CFG["scf_ckpt"]), "gene", device=device)
    cfg = dict(cfg)
    cfg["seq_len"] = SEQ_LEN
    cfg["gene_num"] = SEQ_LEN
    for sub in ("encoder", "decoder"):
        if isinstance(cfg.get(sub), dict):
            cfg[sub]["seq_len"] = SEQ_LEN

    from pretrainmodels import select_model  # noqa: E402

    model = select_model(cfg)
    sd = torch.load(YB / CFG["scf_ckpt"], map_location="cpu", weights_only=False)
    sd = sd["gene"]["state_dict"]
    sd = {k[len("model.") :] if k.startswith("model.") else k: v for k, v in sd.items()}
    sd.pop("pos_emb.weight", None)  # 酵母词表形状不同,加载后整体替换
    missing, unexpected = model.load_state_dict(sd, strict=False)
    bad = [k for k in missing if "pos_emb" not in k]
    assert (
        not bad and not unexpected
    ), f"权重对齐失败 missing={bad[:5]} unexpected={unexpected[:5]}"

    if route in ("A2", "C2"):
        init = np.load(ASSETS / f"{route}_init.npy")
        assert init.shape == (SEQ_LEN + 1, 768), init.shape
        emb = nn.Embedding(SEQ_LEN + 1, 768)
        with torch.no_grad():
            emb.weight.copy_(torch.from_numpy(init))
        model.pos_emb = emb
    else:
        prot = np.load(ASSETS / "B2_esm2_matrix.npy")
        model.pos_emb = ProteinEmbeddingInjector(prot, n_special=3, d_model=768)

    enc_layers = model.encoder.transformer_encoder
    enc_layers = enc_layers.layers if hasattr(enc_layers, "layers") else enc_layers
    for layer in enc_layers[:-2]:
        for p in layer.parameters():
            p.requires_grad = False
    model.encoder = Checkpointed(model.encoder)
    model.decoder = Checkpointed(model.decoder)
    model = model.to(device)
    return model, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "configs/b2_training.json")
    ap.add_argument("--route", default="B2", choices=["B2"])
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument(
        "--smoke-batch",
        type=int,
        default=2,
        help="Number of cells for --smoke only; full training keeps frozen batch 32",
    )
    ap.add_argument("--device", default="auto")
    ap.add_argument(
        "--smoke", action="store_true", help="single train step + table save, then exit"
    )
    args = ap.parse_args()

    outdir = args.output.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        ap.error("Output must be a new or empty directory")
    if args.smoke_batch < 1:
        ap.error("--smoke-batch must be positive")
    outdir.mkdir(parents=True, exist_ok=True)
    from project_io import sha256, write_json

    for key, item in CFG["required_assets"].items():
        path = ROOT / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            ap.error(f"Missing or changed training asset: {key}: {path}")
    write_json(outdir / "run_config.json", CFG)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device = pick_device(args.device)
    print(f"[setup] route={args.route} device={device}", flush=True)

    counts = np.load(ASSETS / "corpus_counts.npy", mmap_mode="r")
    n_cells = counts.shape[0]
    model, cfg = build_model(args.route, device)
    trainable = [p for p in model.parameters() if p.requires_grad]
    print(
        f"[setup] cells={n_cells} trainable_params={sum(p.numel() for p in trainable)/1e6:.1f}M",
        flush=True,
    )
    opt = torch.optim.AdamW(trainable, lr=LR)

    rng = np.random.default_rng(SEED)
    perm = rng.permutation(n_cells)
    n_val = int(0.05 * n_cells)
    val_idx, train_idx = perm[:n_val], perm[n_val:]
    import csv

    with (ASSETS / "corpus_cells.tsv").open() as f:
        cells = list(csv.DictReader(f, delimiter="\t"))
    assert len(cells) == n_cells
    validation = set(val_idx.tolist())
    with (outdir / "split.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_index", "cell", "split"])
        for i, cell in enumerate(cells):
            writer.writerow(
                [i, cell["cell"], "validation" if i in validation else "train"]
            )

    def build_subset(batch_raw):
        """每细胞: 表达基因中均匀随机抽 MAX_SEQ 个(trunc_by_sample 官方行为),
        附 2 个分辨率位; 返回 (x_sub, gene_ids_sub), 位置 id = 真实基因列号。"""
        raw = torch.from_numpy(np.asarray(batch_raw, dtype=np.float32))
        total = raw.sum(dim=1, keepdim=True).clamp(min=1.0)
        vals = torch.log1p(raw / total * 1e4)
        res_v = torch.cat(
            [torch.full((raw.shape[0], 1), 4.0), torch.log10(total)], dim=1
        )
        xs, gs = [], []
        for b in range(vals.shape[0]):
            idx = torch.nonzero(vals[b] > 0, as_tuple=True)[0]
            if idx.numel() > MAX_SEQ:
                perm = torch.randperm(idx.numel())[:MAX_SEQ]
                idx = idx[perm]
            xs.append(torch.cat([vals[b, idx], res_v[b]]))
            gs.append(torch.cat([idx, torch.tensor([N_GENES, N_GENES + 1])]))
        K = max(x.numel() for x in xs)
        x_sub = torch.zeros((len(xs), K), dtype=torch.float32)
        g_sub = torch.full(
            (len(gs), K), SEQ_LEN, dtype=torch.long
        )  # 尾部填充位 → pad 行
        for i, (xv, gv) in enumerate(zip(xs, gs)):
            x_sub[i, : xv.numel()] = xv
            g_sub[i, : gv.numel()] = gv
        return x_sub, g_sub

    def step(batch_raw, train):
        x, gene_ids = build_subset(batch_raw)
        x, gene_ids = x.to(device), gene_ids.to(device)
        if FIXED:
            # Actual gene IDs identify eligible positions independently in each row.
            eligible = gene_ids < N_GENES
            r = torch.rand_like(x)
            mask = eligible & (
                ((x > 0) & (r < MASK_P)) | ((x == 0) & (r < ZERO_MASK_P))
            )
            if not mask.any():
                mask.flatten()[
                    torch.nonzero(eligible.flatten(), as_tuple=True)[0][0]
                ] = True
            masked = x.clone()
            masked[mask] = float(cfg["mask_token_id"])
        else:
            expressed = x[:, :-2] > 0
            r = torch.rand_like(x[:, :-2])
            mask = (expressed & (r < MASK_P)) | (~expressed & (r < ZERO_MASK_P))
            masked = x.clone()
            masked[:, :-2][mask] = float(cfg["mask_token_id"])
        enc = encdec_subset(masked, x, gene_ids, cfg)
        (enc_data, enc_pos, enc_pad, enc_labels, dec_data, dec_pad, _, _, dec_pos) = enc
        with torch.autocast(
            "cuda", dtype=torch.bfloat16, enabled=device.type == "cuda"
        ):
            out = model(
                x=enc_data,
                padding_label=enc_pad,
                encoder_position_gene_ids=enc_pos,
                encoder_labels=enc_labels,
                decoder_data=dec_data,
                mask_gene_name=False,
                mask_labels=None,
                decoder_position_gene_ids=dec_pos,
                decoder_data_padding_labels=dec_pad,
            )
            loss = (
                nn.functional.mse_loss(out[mask], x[mask])
                if FIXED
                else nn.functional.mse_loss(out[:, :-2][mask], x[:, :-2][mask])
            )
        if train:
            (loss / ACCUM).backward()
        return float(loss.detach()), int(mask.sum())

    def save_table(tag):
        model.eval()
        with torch.no_grad():
            if args.route == "B2":
                inj = model.pos_emb
                table = inj.proj(inj.prot).cpu().numpy().astype(np.float32)
            else:
                table = (
                    model.pos_emb.weight[:N_GENES]
                    .detach()
                    .cpu()
                    .numpy()
                    .astype(np.float32)
                )
        np.save(outdir / f"gene_table_{tag}.npy", table)
        model.train()

    t0 = time.time()
    if args.smoke:
        batch = counts[np.sort(train_idx[: args.smoke_batch])]
        loss, n_masked = step(batch, train=True)
        nn.utils.clip_grad_norm_(trainable, CLIP)
        opt.step()
        save_table("smoke")
        print(f"[smoke] one step OK, loss={loss:.4f}; table saved", flush=True)
        return
    metrics = []
    for epoch in range(1, EPOCHS + 1):
        epoch_started = time.time()
        model.train()
        order = rng.permutation(len(train_idx))
        opt.zero_grad(set_to_none=True)
        running, seen = 0.0, 0
        train_sse, train_n = 0.0, 0
        for b in range(0, len(order), BATCH):
            batch = counts[np.sort(train_idx[order[b : b + BATCH]])]
            loss, n_masked = step(batch, train=True)
            running += loss
            train_sse += loss * n_masked
            train_n += n_masked
            seen += 1
            if seen % ACCUM == 0:
                nn.utils.clip_grad_norm_(trainable, CLIP)
                opt.step()
                opt.zero_grad(set_to_none=True)
            if seen % 200 == 0:
                el = time.time() - t0
                print(
                    f"[ep{epoch}] step {seen}/{len(order)//BATCH} loss={running/seen:.4f} "
                    f"elapsed={el/60:.1f}min",
                    flush=True,
                )
        vloss, vb, val_sse, val_n = 0.0, 0, 0.0, 0
        if FIXED:
            model.eval()
        # Fork preserves the training RNG stream while replaying exactly the same
        # gene subsampling and validation masks at every epoch.
        import contextlib

        with (
            torch.random.fork_rng(
                devices=[device.index or 0] if device.type == "cuda" else []
            )
            if FIXED
            else contextlib.nullcontext()
        ):
            if FIXED:
                torch.manual_seed(SEED + 100000)
            with torch.no_grad() if FIXED else contextlib.nullcontext():
                for b in range(0, len(val_idx), BATCH):
                    loss, n_masked = step(
                        counts[np.sort(val_idx[b : b + BATCH])], train=False
                    )
                    vloss += loss
                    vb += 1
                    val_sse += loss * n_masked
                    val_n += n_masked
        metrics.append(
            {
                "epoch": epoch,
                "train_mse_batch_mean": running / seen,
                "validation_mse_batch_mean": vloss / vb,
                "train_mse_mask_weighted": train_sse / train_n,
                "validation_mse_mask_weighted": val_sse / val_n,
                "train_masked_entries": train_n,
                "validation_masked_entries": val_n,
                "seconds": time.time() - epoch_started,
            }
        )
        write_json(outdir / "epoch_metrics.json", metrics)
        print(
            f"[ep{epoch}] train={running/max(seen,1):.4f} val={vloss/max(vb,1):.4f}",
            flush=True,
        )
        save_table(f"ep{epoch}")

    torch.save(
        {
            "state_dict": {k: v.cpu() for k, v in model.state_dict().items()},
            "route": args.route,
            "implementation": CFG.get("implementation", "historical-v1"),
            "training_config": CFG,
            "epoch_metrics": metrics,
            "config": cfg,
        },
        outdir / "final_model.pt",
    )
    save_table("final")
    meta = {
        "route": args.route,
        "epochs": EPOCHS,
        "batch": BATCH,
        "accum": ACCUM,
        "lr": LR,
        "mask_p": MASK_P,
        "zero_mask_p": ZERO_MASK_P,
        "seed": SEED,
        "n_train": len(train_idx),
        "n_val": len(val_idx),
        "implementation": CFG.get("implementation", "historical-v1"),
        "validation_fixed_masks": FIXED,
        "checkpoint_sha256": sha256(outdir / "final_model.pt"),
        "gene_table_sha256": sha256(outdir / "gene_table_final.npy"),
        "split_sha256": sha256(outdir / "split.tsv"),
    }
    (outdir / "train_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[done] {outdir} ({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
