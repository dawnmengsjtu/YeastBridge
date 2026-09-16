#!/usr/bin/env python3
"""ej_dc_nominate — EJ-dc 两阶段确认的机械提名（规则冻结于 ESM2_JOINT_DC_ADDENDUM.md §3）。

输入：EJ-dc 双向全屏结果目录（exec_matrix.tsv）。
规则（冻结）：每方向按 z=(rho-null_mean)/null_sd 降序 top40；跨方向按对去重；
输出 pair_allowlist_{pos,neg}40.tsv + nomination_record.json（含阶段一 sha、
z 范围、化合物集中度、家族结构注记）。
用法：
  python scripts/ej_dc_nominate.py --pos-screen <dir> --neg-screen <dir> --out <panel_dir>
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pos-screen", required=True)
    ap.add_argument("--neg-screen", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=40)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rec = {"mode": "ej_dc_two_stage_nomination",
           "rule_doc": "docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md (section 3)",
           "rule": f"top-{args.top_n} by z per direction, dedup across directions",
           "stage1_screens": {}, "directions": {}}
    frames = {}
    for tag, d in (("pos", args.pos_screen), ("neg", args.neg_screen)):
        p = Path(d) / "exec_matrix.tsv"
        rec["stage1_screens"][tag] = {"path": str(p), "exec_matrix_sha256": sha(p)}
        df = pd.read_csv(p, sep="\t",
                         usecols=["target_id", "inchikey", "spearman_rho",
                                  "emp_p", "null_mean", "null_sd"])
        df["z"] = (df.spearman_rho - df.null_mean) / df.null_sd.replace(0, np.nan)
        top = df.nlargest(args.top_n, "z")
        top[["target_id", "inchikey"]].to_csv(out / f"pair_allowlist_{tag}40.tsv",
                                              sep="\t", index=False)
        frames[tag] = top
        n_floor = int((df.emp_p <= 1.0 / 1001 + 1e-12).sum())
        rec["stage1_screens"][tag]["floor_pairs"] = n_floor
        rec["stage1_screens"][tag]["floor_null_expectation"] = round(len(df) / 1001, 1)
        rec["directions"][tag] = {
            "n": len(top), "z_range": [float(top.z.min()), float(top.z.max())],
            "unique_compounds": int(top.inchikey.nunique()),
            "compound_counts_top4": top.inchikey.value_counts().head(4).to_dict()}
        print(tag, "top written; compounds:", top.inchikey.nunique(),
              "floor:", n_floor, "vs expect", round(len(df) / 1001, 1))
    pos, neg = frames["pos"][["target_id", "inchikey"]], frames["neg"][["target_id", "inchikey"]]
    rec["cross_direction_overlap_pairs"] = int(len(pd.merge(pos, neg, on=["target_id", "inchikey"])))
    for f in ("pair_allowlist_pos40.tsv", "pair_allowlist_neg40.tsv"):
        rec[f + ".sha256"] = sha(out / f)
    (out / "nomination_record.json").write_text(json.dumps(rec, indent=1))
    print("[done] ->", out)


if __name__ == "__main__":
    main()
