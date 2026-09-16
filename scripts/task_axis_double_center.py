#!/usr/bin/env python3
"""task_axis_double_center — EJ 任务轴共模去除（EJ-dc，2026-09-12 用户指示优化）。

动机（数据在案，results_esm2joint_screen_20260912）：+z 全屏 top25 提名全部是
同一化合物（SQFNCEKTEKNPCN, dose 12.2）× 25 个互不相关靶点、rho≈0.052 同值；
emp_p 达底对 2100 < 全局零假设期望 ~3800 —— 任务排序共享菌株共模，被泛化
化合物谱批量假命中，尾部无真信号。共模同时抬高假阳性并吞掉真对的排序方差。

变换（确定性后处理 EJ 冻结导出分，不重训任何算子）：
  1. S_pool = 各靶点 ensemble_score 在 1279 池株上的矩阵 (targets x strains)
  2. S1 = S_pool - 列均值（每株跨全部靶点去均值）      —— 去菌株共模
  3. S2 = S1 - PC1（S1 的 SVD 首分量，双侧投影去除）    —— 去残余共享主轴
  4. 每靶点在 S2 上重排；尾段 3 株按 ORF 字典序单调拼接（与 EJ 导出同规）
记录：共模方差占比（列均值 + PC1）、输入目录逐文件哈希摘要、输出哈希。

用法：
  python scripts/task_axis_double_center.py \
      --input <stage>/esm2_joint_tasks_20260912 \
      --output <stage>/esm2_joint_tasks_dc_20260912
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    indir, outdir = Path(args.input), Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    rec_in = json.loads((indir / "export_record.json").read_text())
    pool_n = rec_in["pool"]
    tail = list(rec_in["tail_strains"])

    files = sorted(indir.glob("yeast_task_*.tsv"))
    print(f"task files: {len(files)}", flush=True)
    targets, rows = [], []
    orfs = None
    for fp in files:
        df = pd.read_csv(fp, sep="\t")
        sc = dict(zip(df.yeast_gene, df.ensemble_score))
        if orfs is None:
            orfs = list(df.yeast_gene)  # 基准株集合(首文件 rank 序, 仅作行索引)
        assert set(sc) == set(orfs), f"strain set mismatch: {fp}"
        rows.append([sc[o] for o in orfs])
        targets.append(fp.name[len("yeast_task_"):-len(".tsv")])
    S_full = np.asarray(rows, dtype=np.float64)          # (T, 1282)
    pool_mask = np.array([o not in set(tail) for o in orfs])
    assert pool_mask.sum() == pool_n, (pool_mask.sum(), pool_n)
    S = S_full[:, pool_mask]                              # (T, 1279)

    total_var = float((S ** 2).sum())
    colmean = S.mean(axis=0, keepdims=True)
    var_colmean = float((colmean ** 2).sum() * S.shape[0]) / total_var
    S1 = S - colmean
    U, sv, Vt = np.linalg.svd(S1, full_matrices=False)
    pc1_var = float(sv[0] ** 2) / float((sv ** 2).sum())
    S2 = S1 - np.outer(U[:, 0] * sv[0], Vt[0])
    print(f"common-mode variance: strain-mean {var_colmean:.3f}, "
          f"residual PC1 {pc1_var:.3f}", flush=True)

    # 每靶点重排 + 尾段拼接（与 EJ 导出同规）
    report = {"mode": "task_axis_double_center",
              "input_dir": str(indir),
              "input_export_record_sha256": sha256_file(str(indir / "export_record.json")),
              "transform": ["colmean_across_targets", "pc1_removal"],
              "common_mode_variance": {"strain_mean_frac": var_colmean,
                                       "residual_pc1_frac": pc1_var},
              "n_targets": len(targets), "pool": pool_n, "tail_strains": tail,
              "targets": {}}
    tail_idx = [i for i, o in enumerate(orfs) if o in set(tail)]
    tail_order = sorted(range(len(tail_idx)), key=lambda a: orfs[tail_idx[a]])
    for ti, t in enumerate(targets):
        s = S2[ti]
        full = np.empty(S_full.shape[1])
        full[pool_mask] = s
        if tail_idx:
            lo = float(s.min())
            for rank_a, idx_a in enumerate(tail_order):
                full[tail_idx[idx_a]] = lo - 0.001 - rank_a * 1e-9
        order = np.argsort(-full, kind="stable")
        fp = outdir / f"yeast_task_{t}.tsv"
        with fp.open("w", encoding="utf-8", newline="") as fh:
            fh.write("rank\tyeast_gene\tensemble_score\n")
            for rk, j in enumerate(order, 1):
                fh.write(f"{rk}\t{orfs[j]}\t{full[j]:.9f}\n")
        report["targets"][t] = {"status": "OK", "rows": len(orfs), "path": str(fp)}
        if (ti + 1) % 200 == 0:
            print(f"  rewritten {ti + 1}/{len(targets)}", flush=True)
    (outdir / "export_record.json").write_text(json.dumps(report, indent=1) + "\n")
    print(f"[done] -> {outdir}", flush=True)


if __name__ == "__main__":
    main()
