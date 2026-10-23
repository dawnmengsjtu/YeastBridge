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


def project_genetic_interactions(scores, genes, archive):
    """Transfer correspondence weights through measured SGA effects, with missingness explicit."""
    import zipfile
    genes = np.asarray(genes)
    position = {g: i for i, g in enumerate(genes)}
    total = np.zeros((len(genes), len(genes)), dtype=np.float64)
    count = np.zeros(total.shape, dtype=np.uint8)
    records = []
    with zipfile.ZipFile(archive) as source:
        for member in sorted(n for n in source.namelist() if n.endswith('.cdt')):
            with source.open(member) as handle:
                header = [handle.readline().decode().rstrip('\r\n').split('\t') for _ in range(6)]
                frame = pd.read_csv(handle, sep='\t', header=None)
            values = frame.iloc[:, 6:].astype(float)
            values.index = frame.iloc[:, 2].astype(str)
            values.columns = header[2][6:]
            values = values.loc[values.index.isin(position), values.columns.isin(position)]
            values = values.groupby(level=0).mean().T.groupby(level=0).mean().T
            matrix = values.to_numpy(); valid = np.isfinite(matrix)
            ix = np.ix_([position[g] for g in values.index], [position[g] for g in values.columns])
            total[ix] += np.where(valid, matrix, 0.0)
            count[ix] += valid.astype(np.uint8)
            records.append(dict(member=member, source_gene_rows=len(values), source_gene_columns=len(values.columns), observed_cells=int(valid.sum())))
    measured = count > 0
    directed = np.divide(total, count, out=np.zeros_like(total), where=measured)
    orientations = measured.astype(np.uint8) + measured.T.astype(np.uint8)
    effects = np.divide(directed + directed.T, orientations, out=np.zeros_like(total), where=orientations > 0)
    measured = orientations > 0
    np.fill_diagonal(measured, False)
    effects[~measured] = 0
    if np.any(scores < 0) or not np.isfinite(scores).all():
        raise ValueError('genetic projection requires uncorrected nonnegative correspondence weights')
    weights = scores / scores.sum(axis=1, keepdims=True)
    denominator = weights @ measured.astype(np.float64)
    keep = np.all(denominator > 0, axis=0)
    projected = -(weights @ effects)[:, keep] / denominator[:, keep]
    if not np.isfinite(projected).all() or int(keep.sum()) < 100:
        raise ValueError('invalid genetic-interaction projection')
    report = dict(archive=str(archive), sha256=sha256_file(archive), source_doi='10.1126/science.aaf1420',
                  source_blocks=records, input_strains=len(genes), output_strains=int(keep.sum()),
                  unobserved_strains=genes[~keep].tolist(), observed_undirected_pairs=int(np.triu(measured, 1).sum()),
                  projection='negative weighted SGA epsilon divided by available correspondence-weight mass; no human bioactivity labels',
                  alleles='mean across row alleles, then column alleles; source blocks and measured orientations averaged; self-interactions excluded')
    return projected, genes[keep].tolist(), report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--remove-pcs", type=int, choices=(0, 1), default=1,
                    help="Frozen sensitivity comparison: strain-mean only or strain-mean plus PC1")
    ap.add_argument("--genetic-interaction-archive")
    ap.add_argument("--projection-raw-output")
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

    projection_record = None
    if args.genetic_interaction_archive:
        if not args.projection_raw_output or tail:
            ap.error('genetic projection requires a raw-output directory and an input without an artificial tail')
        S, orfs, projection_record = project_genetic_interactions(S, orfs, args.genetic_interaction_archive)
        S_full = S; pool_n = len(orfs); pool_mask = np.ones(pool_n, dtype=bool)
        rawdir = Path(args.projection_raw_output); rawdir.mkdir(parents=True, exist_ok=True)
        rawrec = dict(mode='genetic_interaction_projection', input_dir=str(indir), input_export_record_sha256=sha256_file(indir/'export_record.json'),
                      n_targets=len(targets), pool=pool_n, tail_strains=[], genetic_interaction_projection=projection_record, targets={})
        for ti, target in enumerate(targets):
            order = np.argsort(-S[ti], kind='stable'); fp = rawdir/f'yeast_task_{target}.tsv'
            pd.DataFrame(dict(rank=np.arange(1, pool_n+1), yeast_gene=np.asarray(orfs)[order], ensemble_score=S[ti, order])).to_csv(fp, sep='\t', index=False, float_format='%.9f')
            rawrec['targets'][target] = dict(status='OK', rows=pool_n, path=str(fp))
        (rawdir/'export_record.json').write_text(json.dumps(rawrec, indent=1))
        print('genetic-interaction projection', json.dumps(projection_record), flush=True)

    total_var = float((S ** 2).sum())
    colmean = S.mean(axis=0, keepdims=True)
    var_colmean = float((colmean ** 2).sum() * S.shape[0]) / total_var
    S1 = S - colmean
    U, sv, Vt = np.linalg.svd(S1, full_matrices=False)
    pc1_var = float(sv[0] ** 2) / float((sv ** 2).sum())
    S2 = S1 - np.outer(U[:, 0] * sv[0], Vt[0]) if args.remove_pcs else S1
    print(f"common-mode variance: strain-mean {var_colmean:.3f}, "
          f"residual PC1 {pc1_var:.3f}", flush=True)

    # 每靶点重排 + 尾段拼接（与 EJ 导出同规）
    report = {"mode": "task_axis_double_center",
              "input_dir": str(indir),
              "input_export_record_sha256": sha256_file(str(indir / "export_record.json")),
              "transform": ["colmean_across_targets"] + (["pc1_removal"] if args.remove_pcs else []),
              "n_removed_pcs": args.remove_pcs,
              "common_mode_variance": {"strain_mean_frac": var_colmean,
                                       "residual_pc1_frac": pc1_var},
              "n_targets": len(targets), "pool": pool_n, "tail_strains": tail,
              "targets": {}}
    if projection_record is not None:
        report['genetic_interaction_projection'] = projection_record
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
