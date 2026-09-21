#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""predict.py — 一键生成标准化候选清单 results.csv（附件5 第四节）。

输入：本仓库 results_a6000/ 内的两阶段确认与家族特异性在册结果
      + data/compounds_smiles.tsv（80 对涉及化合物的 SMILES/CID 对照）。
输出：results.csv（UTF-8），字段：候选编号/赛道/靶点/化合物/SMILES/
      方向/关键指标(spearman_rho, emp_p, q_bh)/证据等级/剂量/模型版本/备注。
用法：python predict.py  （工作目录为仓库根）
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
EH = ROOT / "results_a6000" / "execute_hiphop"
TRACK = "赛道三（离子通道、GPCR 等药靶的小分子药物虚拟筛选）"
MODEL = "YeastBridge EJ-dc 任务轴（RRF20: direct_pc1 + aligned_C + b2_verbatim，共模双去除）"


def main():
    pos = pd.read_csv(EH / "results_esm2jointdc_confirm100k_pos_20260912/exec_matrix.tsv", sep="\t")
    neg = pd.read_csv(EH / "results_esm2jointdc_neg_confirm100k_neg_20260912/exec_matrix.tsv", sep="\t")
    pos["direction"], neg["direction"] = "+z（耐药富集）", "-z（超敏）"

    # 靶点级证据：残差检验过阈（p<6.25e-4，在册口径）
    resid_sig = set()
    for d in ("results_esm2jointdc_resid_pos_20260913",
              "results_esm2jointdc_resid_neg_neg_20260913"):
        r = pd.read_csv(EH / d / "exec_matrix.tsv", sep="\t")
        resid_sig |= set(map(tuple, r[r.emp_p < 6.25e-4][["target_id", "inchikey"]].values))
    # 构造性靶点级：0.9 阈值下无近同谱的单例靶点（在册 build_record）
    build = json.loads((ROOT / "panels/family_specificity_panel_20260913/"
                        "build_record.json").read_text())
    singles = set(build.get("nominated_singletons", []))

    comp = pd.read_csv(ROOT / "data" / "compounds_smiles.tsv", sep="\t")
    cmap = dict(zip(comp.inchikey, comp.smiles))
    cinfo = comp.set_index("inchikey")[["pubchem_cid"]].to_dict("index")

    rows = pd.concat([pos, neg], ignore_index=True).sort_values(["direction", "target_id"])
    out = []
    for i, r in enumerate(rows.itertuples(), 1):
        if (r.target_id, r.inchikey) in resid_sig:
            tier = "靶点级（残差检验过阈）"
            note = "去掉家族公共分量后仍显著"
        elif r.target_id in singles:
            tier = "靶点级（构造性，无近同谱）"
            note = "0.9 阈值下为单例靶点"
        else:
            tier = "家族级"
            note = "家族均值轴复现（12/12 在册）"
        out.append({
            "候选编号": f"YB-{i:03d}", "赛道": TRACK, "靶点": r.target_id,
            "化合物_InChIKey": r.inchikey,
            "化合物_CID": cinfo.get(r.inchikey, {}).get("pubchem_cid", ""),
            "SMILES": cmap.get(r.inchikey, ""), "方向": r.direction,
            "spearman_rho": round(r.spearman_rho, 4),
            "emp_p": round(r.emp_p, 6), "q_bh": round(r.q, 8),
            "剂量": r.dose, "证据等级": tier, "模型与版本": MODEL,
            "备注": note,
        })
    df = pd.DataFrame(out)
    df.to_csv(ROOT / "results.csv", index=False, encoding="utf-8")
    print(f"results.csv 生成：{len(df)} 条候选（靶点级 "
          f"{(df.证据等级 != '家族级').sum()} / 家族级 {(df.证据等级 == '家族级').sum()}）")


if __name__ == "__main__":
    main()
