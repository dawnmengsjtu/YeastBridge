#!/usr/bin/env python3
"""esm2_joint_task_export — EJ 产品线任务轴导出（V9 方法冻结配方）。

方法冻结：docs/cross_species_match/ESM2_JOINT_METHOD_FREEZE.md
配置：configs/esm2_joint_tasks.json（universe 资产哈希冻结；其余资产经
v7_setup 按 configs/v7_b2_joint_formal.json 哈希核对）

每靶点任务分 = RRF k=20 三成员：
  m1 direct_pc1     cos(U1, Y1)          [V8 臂 B 注册 0.9638]
  m2 aligned_C      npc*=3 空间三算子秩平均（ridge λ1 + anchored λ1 + PLS8,
                    官方 train 对训练）   [V8 臂 C 注册 0.8770, 门 A 过]
  m3 b2_verbatim    cos(l2n(U@W_inj+b), l2n(routeB))  [B2 原版, 跨基座成员]
超参运行时断言与 V8 冻结选择一致。池 1282 中 3 株缺表示 → ORF 字典序尾段。

用法：
  python scripts/esm2_joint_task_export.py --config configs/esm2_joint_tasks.json \
      --output <stage>/esm2_joint_tasks_20260912
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v5_joint_formal_run import l2n, sha256_file  # noqa: E402
from v7_b2_formal_run import v7_setup  # noqa: E402
from v8_esm2_formal_run import apply_pc, fit_pc  # noqa: E402


def rank_desc(v):
    order = np.argsort(-v, kind="stable")
    rk = np.empty(len(v), dtype=np.int64)
    rk[order] = np.arange(1, len(v) + 1)
    return rk


def pct_rank(v):
    return (np.argsort(np.argsort(v)) + 1) / len(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    root = Path(cfg["stage_root"])
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    mp = cfg["method"]

    # ---- universe 资产哈希核对 ----
    manifest = {}
    for k, v in cfg["inputs"].items():
        if not isinstance(v, dict) or "path" not in v or "sha256" not in v:
            continue
        p = root / v["path"]
        h = sha256_file(str(p))
        manifest[k] = {"path": v["path"], "sha256": h, "match": h == v["sha256"]}
        if h != v["sha256"]:
            raise SystemExit(f"ABORT: hash mismatch {k}")
        if "index_path" in v:
            hi = sha256_file(str(root / v["index_path"]))
            manifest[k + ".index"] = {"path": v["index_path"], "sha256": hi,
                                      "match": hi == v["index_sha256"]}
            if hi != v["index_sha256"]:
                raise SystemExit(f"ABORT: index hash mismatch {k}")
    print("universe asset hashes: all match", flush=True)

    # ---- V7 机械复用（池/配对/W_inj/酵母表, 哈希核对在内） ----
    # v7_config 路径相对仓库根（yeastbridge_re_mvp），非 stage_root
    st = v7_setup(str(HERE.parent / mp["v7_config"]))

    # ---- 超参断言：与 V8 冻结选择一致 ----
    v8res = json.loads((root / cfg["inputs"]["v8_result_for_hyper_assert"]["path"]).read_text())
    info = v8res["trackA"]["armC_recipe_info"]
    assert info["npc_selected"] == mp["npc_aligned"], (info["npc_selected"], mp["npc_aligned"])
    assert float(info["picks"]["ridge"]["hyper"]) == float(mp["ridge_lambda"])
    assert float(info["picks"]["anchored"]["hyper"]) == float(mp["anchored_lambda"])
    assert int(info["picks"]["pls"]["hyper"]) == int(mp["pls_k"])
    print("hyper assert vs V8 frozen picks: OK "
          f"(npc={mp['npc_aligned']}, lam_r={mp['ridge_lambda']}, "
          f"lam_a={mp['anchored_lambda']}, pls_k={mp['pls_k']})", flush=True)

    # ---- universe 嵌入 ----
    U = np.load(root / cfg["inputs"]["universe_esm2"]["path"]).astype(np.float64)
    uidx = pd.read_csv(root / cfg["inputs"]["universe_esm2"]["index_path"],
                       sep="\t").fillna("")
    symbols = [str(s) for s in uidx["common"]]
    assert U.shape[0] == len(symbols)
    seen, rows_keep = set(), []
    for i, s in enumerate(symbols):
        if s and s not in seen:
            seen.add(s)
            rows_keep.append(i)
    n_dup = len(symbols) - len(rows_keep)
    symbols = [symbols[i] for i in rows_keep]
    U = l2n(U[rows_keep])
    print(f"universe targets: {len(symbols)} (dup/blank dropped {n_dup})", flush=True)

    # ---- 冻结 PC 基（与 V8/V9 同拟合矩阵） ----
    fx = fit_pc(st["X_es_all"])
    fy = fit_pc(st["Y_es"])
    U1 = apply_pc(U, fx[0], fx[1], mp["npc_direct"])
    Y1 = apply_pc(st["Y_es"], fy[0], fy[1], mp["npc_direct"])
    UN = apply_pc(U, fx[0], fx[1], mp["npc_aligned"])
    YN = apply_pc(st["Y_es"], fy[0], fy[1], mp["npc_aligned"])

    # ---- m2 算子（官方 train 对, 冻结超参, 确定性重拟合） ----
    prim, q_pos, pool_pos = st["prim"], st["q_pos"], st["pool_pos"]
    q_rows = np.array([q_pos[h] for h, _, _, _ in prim])
    y_rows = np.array([pool_pos[y] for _, y, _, _ in prim])
    tr = np.array([s == "train" for _, _, s, _ in prim])
    hrow = [st["sym2row_h"][g] for g in st["query_genes"]]
    XN = apply_pc(st["X_es_all"][hrow], fx[0], fx[1], mp["npc_aligned"])
    Xtr, Ytr = XN[q_rows[tr]], YN[y_rows[tr]]
    XtX, XtY = Xtr.T @ Xtr, Xtr.T @ Ytr
    I = np.eye(XN.shape[1])
    W_r = np.linalg.solve(XtX + mp["ridge_lambda"] * I, XtY)
    W_a = np.linalg.solve(XtX + mp["anchored_lambda"] * I,
                          XtY + mp["anchored_lambda"] * I)
    mx_, sx_ = Xtr.mean(0), np.where(Xtr.std(0) > 0, Xtr.std(0), 1)
    my_, sy_ = Ytr.mean(0), np.where(Ytr.std(0) > 0, Ytr.std(0), 1)
    Usv, _, Vsv = np.linalg.svd(((Xtr - mx_) / sx_).T @ ((Ytr - my_) / sy_),
                                full_matrices=False)
    k = int(mp["pls_k"])
    Px_u = l2n(((UN - mx_) / sx_) @ Usv[:, :k])
    Py = l2n(((YN - my_) / sy_) @ Vsv[:k].T)
    print(f"operators fitted on {int(tr.sum())} train pairs", flush=True)

    # ---- m3 B2 verbatim ----
    Q_b2 = l2n(U @ st["W_inj"].T + st["b_inj"])
    Yrb = st["Y_rb"]

    # ---- 池与尾段 ----
    orfs_all = [str(x) for x in st["meta_y"].target_stable_id]
    pool = st["pool"]
    orf_row = {o: i for i, o in enumerate(orfs_all)}
    missing = sorted(o for o in orfs_all if o not in pool_pos)
    print(f"pool {len(pool)}/{len(orfs_all)}, tail strains: {missing}", flush=True)

    kR = mp["rrf_k"]
    report = {"mode": "export_esm2_joint_tasks",
              "method_freeze": "ESM2_JOINT_METHOD_FREEZE.md",
              "config_sha256": sha256_file(args.config),
              "script_sha256": sha256_file(str(Path(__file__).resolve())),
              "members": mp["members"], "rrf_k": kR,
              "hyper_frozen": {"npc_direct": mp["npc_direct"],
                               "npc_aligned": mp["npc_aligned"],
                               "ridge_lambda": mp["ridge_lambda"],
                               "anchored_lambda": mp["anchored_lambda"],
                               "pls_k": mp["pls_k"]},
              "train_pairs": int(tr.sum()), "pool": len(pool),
              "tail_strains": missing, "n_targets": len(symbols),
              "targets": {}}
    for i, t in enumerate(symbols):
        s1 = Y1 @ U1[i]
        Sr = YN @ (UN[i] @ W_r)
        Sa = YN @ (UN[i] @ W_a)
        Sp = Py @ Px_u[i]
        s2 = (pct_rank(Sr) + pct_rank(Sa) + pct_rank(Sp)) / 3.0
        s3 = Yrb @ Q_b2[i]
        rrf = (1.0 / (kR + rank_desc(s1)) + 1.0 / (kR + rank_desc(s2))
               + 1.0 / (kR + rank_desc(s3)))
        full = np.empty(len(orfs_all))
        for o, j in pool_pos.items():
            full[orf_row[o]] = rrf[j]
        if missing:
            lo = float(rrf.min())
            for a_, o in enumerate(missing):
                full[orf_row[o]] = lo - 0.001 - a_ * 1e-9
        order = np.argsort(-full, kind="stable")
        fp = outdir / f"yeast_task_{t}.tsv"
        with fp.open("w", encoding="utf-8", newline="") as fh:
            fh.write("rank\tyeast_gene\tensemble_score\n")
            for rk, j in enumerate(order, 1):
                fh.write(f"{rk}\t{orfs_all[j]}\t{full[j]:.9f}\n")
        report["targets"][t] = {"status": "OK", "rows": len(orfs_all), "path": str(fp)}
        if (i + 1) % 200 == 0:
            print(f"  exported {i + 1}/{len(symbols)}", flush=True)
    (outdir / "export_record.json").write_text(json.dumps(report, indent=1) + "\n")
    print(f"[done] {len(symbols)} targets -> {outdir}", flush=True)


if __name__ == "__main__":
    main()
