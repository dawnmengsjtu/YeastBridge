#!/usr/bin/env python3
"""ESM2 跨物种蛋白匹配执行器（Phase 2c 的工作脚本；编排器只负责调用）。

对任何人源靶点：获取 ESM2 嵌入（universe 预计算 -> 缓存 -> UniProt 取序列 +
extract_esm2_candidates.py 推理），然后给出两臂匹配与代理扰动证据：

- 臂 1 esm2_mean（参照特征）：裸余弦 vs 酵母 ESM2 表；
- 臂 2 B2_scf_esm2inject（scF 在评分链）：ESM2 -> 训练好的 route-B 注入投影
  -> 余弦 vs 训练后 route-B 酵母基因表（机制与 scripts/product_transfer_route_b.py
  逐字一致）；
- 代理证据（现有酵母扰动数据）：每个 top 匹配查 GSE42528 KO 签名 top 响应基因
  与 Lee 2014 HIP/HOP 最优 FD 排名。

参数见 configs/esm2_cross_match.json；输出 JSON 证据卡片段。
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def http_json(url, retries=3):
    import time
    for t in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read())
        except Exception:
            time.sleep(1 + t)
    return None


def human_embedding(target, cfg, py):
    """universe 预计算 -> 缓存 -> UniProt + 推理。返回 (emb, source)。"""
    import numpy as np
    import pandas as pd
    udir = Path(cfg["universe_esm2_dir"])
    if (udir / "index.tsv").exists():
        uidx = pd.read_csv(udir / "index.tsv", sep="\t").fillna("")
        hit = uidx[uidx.common == target]
        if len(hit):
            return np.load(udir / "esm2_mean_fp32.npy")[hit.index[0]], "universe_precomputed"
    cache = Path(cfg["query_cache_dir"]) / target
    if (cache / "esm2_mean_fp32.npy").exists():
        return np.load(cache / "esm2_mean_fp32.npy")[0], "cache"
    q = http_json("https://rest.uniprot.org/uniprotkb/search?query="
                  f"gene_exact:{target}+AND+organism_id:9606+AND+reviewed:true"
                  "&fields=accession,gene_primary,sequence&format=json")
    seq = None
    if q and q.get("results"):
        seq = q["results"][0].get("sequence", {}).get("value")
    if not seq:
        return None, "no_sequence"
    cache.mkdir(parents=True, exist_ok=True)
    acc = q["results"][0].get("primaryAccession", target)
    (cache / "query.fasta").write_text(f">sp|{acc}|{target}_HUMAN GN={target}\n{seq}\n")
    r = subprocess.run(
        [py, str(ROOT / "scripts" / "extract_esm2_candidates.py"),
         "--fasta", str(cache / "query.fasta"), "--outdir", str(cache),
         "--model", cfg["esm2_model"]],
        capture_output=True, text=True, timeout=600)
    if r.returncode != 0 or not (cache / "esm2_mean_fp32.npy").exists():
        return None, f"inference_failed:{r.stderr[-80:]}"
    return np.load(cache / "esm2_mean_fp32.npy")[0], "inferred"


def top_matches(emb, cfg, n):
    """臂 1：裸余弦 vs 酵母 ESM2 表（esm2_mean 参照臂）。"""
    import numpy as np
    import pandas as pd
    idx = pd.read_csv(Path(cfg["yeast_esm2_dir"]) / "index.tsv", sep="\t")
    Y = np.load(Path(cfg["yeast_esm2_dir"]) / "esm2_mean_fp32.npy")
    a = emb / (np.linalg.norm(emb) + 1e-9)
    b = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-9)
    sims = b @ a
    top = sims.argsort()[::-1][:n]
    return [{"systematic": idx.iloc[i]["systematic"], "common": idx.iloc[i]["common"],
             "cosine": round(float(sims[i]), 4)} for i in top]


def route_b_top(emb, cfg, n):
    """臂 2：route B verbatim（scF-B2 注入投影，与 product_transfer_route_b.py 同式）。"""
    import numpy as np
    import pandas as pd
    import torch
    table = np.load(ROOT / cfg["route_b_table"]).astype(np.float64)
    order = pd.read_csv(ROOT / cfg["route_b_gene_order"], sep="\t", dtype=str)["systematic"].tolist()
    ck = torch.load(ROOT / cfg["route_b_model"], map_location="cpu", weights_only=False)
    sd = ck["state_dict"]
    W = sd["pos_emb.proj.weight"].float().numpy()
    b = sd["pos_emb.proj.bias"].float().numpy()
    tn = table / np.maximum(np.linalg.norm(table, axis=1, keepdims=True), 1e-12)
    q = emb.astype(np.float64) @ W.T + b
    q = q / max(np.linalg.norm(q), 1e-12)
    scores = tn @ q
    top = scores.argsort()[::-1][:n]
    idx = pd.read_csv(Path(cfg["yeast_esm2_dir"]) / "index.tsv", sep="\t")
    common_of = dict(zip(idx["systematic"], idx["common"]))
    return [{"systematic": order[i], "common": common_of.get(order[i], order[i]),
             "cosine": round(float(scores[i]), 4)} for i in top]


def proxy_perturbation(matches, cfg, n_sig=5, n_hip=3):
    """现有酵母扰动数据代理证据：GSE42528 KO 签名 + Lee 2014 HIP/HOP。"""
    import pandas as pd
    ym = pd.read_csv(cfg["yeast_matrix"], sep="\t", index_col=0)
    label_of = {str(i).split(":")[-1]: i for i in ym.index}
    fd_files = [f for f in glob.glob(str(Path(cfg["lee2014_fd_dir"]) / "*.txt"))
                if "headers" not in f and "RECOVERY" not in f]
    ev = {}
    for m in matches:
        ent = {"cosine": m["cosine"]}
        lab = label_of.get(m["common"])
        if lab is not None:
            row = ym.loc[lab]
            top = row.abs().sort_values(ascending=False).head(n_sig)
            ent["gse42528_ko"] = {"strain": lab,
                                  "top_responding": {g: round(float(row[g]), 3) for g in top.index}}
        else:
            ent["gse42528_ko"] = "no_strain"
        hits = []
        for f in fd_files:
            try:
                fd = pd.read_csv(f, sep="\t")
                fh = fd[fd.ORF == m["systematic"]]
                if len(fh):
                    # HIP: 超敏 FD 高为证据; HOP: 靶点必需→耐药 FD 低为证据
                    if ".HIP" in f:
                        rank = int((fd.FD > fh.iloc[0].FD).sum()) + 1
                    else:
                        rank = int((fd.FD < fh.iloc[0].FD).sum()) + 1
                    cond = Path(f).name.replace(".HIP.txt", "").replace(".HOP.txt", "")
                    hits.append({"screen": cond, "mode": "HIP" if ".HIP" in f else "HOP",
                                 "fd": round(float(fh.iloc[0].FD), 2),
                                 "rank": rank, "pool": len(fd)})
            except Exception:
                pass
        hits.sort(key=lambda h: h["rank"] / max(h["pool"], 1))
        ent["lee2014_hiphop"] = hits[:n_hip] if hits else "no_data"
        ev[f"{m['systematic']}/{m['common']}"] = ent
    return ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--config", default="configs/esm2_cross_match.json")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--topn", type=int, default=5)
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    py = sys.executable

    emb, src = human_embedding(args.target, cfg, py)
    if emb is None:
        card = {"status": "DEGRADED", "reason": src}
    else:
        arm1 = top_matches(emb, cfg, args.topn)
        arm2 = route_b_top(emb, cfg, args.topn)
        card = {
            "status": "OK", "embedding_source": src, "dim": int(emb.shape[0]),
            "esm2_mean": {"role": "reference_arm", "top_matches": arm1},
            "route_b_scf": {"role": "selected_arm_scf_in_scoring",
                            "mechanism": "ESM2 -> trained B2 injection -> cosine route-B table",
                            "top_matches": arm2},
            "proxy_evidence": proxy_perturbation(arm1, cfg),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(card, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": card["status"], "source": src}, ensure_ascii=False))


if __name__ == "__main__":
    main()
