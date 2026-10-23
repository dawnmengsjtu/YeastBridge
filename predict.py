#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""predict.py — 生成当前离子通道候选清单 results.csv。

输入：ion_channel_scope_20261012/full_measured_assay 的全库匹配及残差结果
      + raw/tier1_response/compounds.tsv.gz 的完整化合物身份表。
输出：results.csv（UTF-8），保留两个完整检验集合各自的全库校正P值。
全部匹配层关联同时保存在分析目录；显式指定历史目录时仍可导出旧分析。
用法：python predict.py  （工作目录为仓库根）
"""
import json
import os
from pathlib import Path

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
EH = ROOT / "results_a6000" / "execute_hiphop"
TRACK = "赛道三（离子通道、GPCR 等药靶的小分子药物虚拟筛选）"
MODEL = "YeastBridge EJ-dc 任务轴（RRF20: direct_pc1 + aligned_C + b2_verbatim，共模双去除）"


def export_screenwide_evidence(run):
    """Export every screen-wide pass and mark the residual-supported subset without a count quota."""
    tables = {}
    for layer, folder in (("match", run), ("residual", run / "residual_fullscreen")):
        null = np.sort(pd.read_csv(folder / "global_tail/screen_max_null.tsv", sep="\t").both.to_numpy())
        plan = json.loads((folder / "global_tail/analysis_plan.json").read_text())
        if len(null) != plan["n"]:
            raise ValueError(f"{layer} full-screen null is incomplete for its recorded plan")
        frames = []
        for direction in ("pos", "neg"):
            frame = pd.read_csv(folder / f"screen_{direction}/exec_matrix.tsv", sep="\t",
                                usecols=["target_id", "inchikey", "dose", "source_condition_index", "spearman_rho", "null_mean", "null_sd"])
            if frame.duplicated(["target_id", "inchikey"]).any() or not frame.null_sd.gt(0).all():
                raise ValueError("invalid full-screen score matrix")
            score = ((frame.spearman_rho - frame.null_mean) / frame.null_sd).to_numpy()
            frame[f"p_{layer}_screenwide"] = (1 + len(null) - np.searchsorted(null, score, side="left")) / (len(null) + 1)
            frame["direction"] = direction
            frames.append(frame)
        tables[layer] = pd.concat(frames, ignore_index=True)
    keys = ["target_id", "inchikey", "direction"]
    full = tables["match"].loc[tables["match"].p_match_screenwide.le(.05)].copy()
    residual = tables["residual"][keys + ["p_residual_screenwide", "dose", "source_condition_index", "spearman_rho"]]
    full = full.merge(residual, on=keys, how="left", suffixes=("", "_residual"), validate="one_to_one")
    full["both_screenwide_tests_pass"] = full.p_residual_screenwide.le(.05)
    # This additional column discloses the two complete-coverage PC settings considered in this round.
    full["both_pass_after_two_variant_bonferroni"] = full.p_match_screenwide.le(.025) & full.p_residual_screenwide.le(.025)
    full["both_pass_after_eight_variant_bonferroni"] = full.p_match_screenwide.le(.05 / 8) & full.p_residual_screenwide.le(.05 / 8)
    full["both_pass_after_twelve_variant_bonferroni"] = full.p_match_screenwide.le(.05 / 12) & full.p_residual_screenwide.le(.05 / 12)
    non_estimable = json.loads((run / "fixed_residual_tasks/export_record.json").read_text())["non_estimable_targets"]
    full["residual_status"] = np.where(full.target_id.isin(non_estimable), "zero_residual_non_estimable",
                                      np.where(full.p_residual_screenwide.notna(), "tested_nonsingleton", "profile_cluster_singleton"))
    scope = pd.read_csv(run.parent / "target_scope.tsv", sep="\t")
    if set(scope.target_family) != {"ion_channel"} or not set(full.target_id) <= set(scope.target_id):
        raise ValueError("screen-wide export is outside the frozen ion-channel scope")
    chemistry = pd.read_csv(os.environ.get("YB_COMPOUND_TABLE", ROOT / "raw/tier1_response/compounds.tsv.gz"), sep="\t")
    chemistry = chemistry[["inchikey", "pubchem_cid", "smiles"]].drop_duplicates()
    if chemistry.inchikey.duplicated().any():raise ValueError("conflicting compound identities")
    full = full.merge(chemistry, on="inchikey", how="left", validate="many_to_one")
    if full.smiles.isna().any():raise ValueError("missing structure for a screened association")
    full["pubchem_cid"] = pd.to_numeric(full.pubchem_cid, errors="raise").astype("Int64")
    full = full.sort_values(["target_id", "inchikey", "direction"]).reset_index(drop=True)
    full.to_csv(run / "screenwide_associations.tsv", sep="\t", index=False)
    supported = full.loc[full.both_screenwide_tests_pass].copy()
    supported.insert(0, "candidate_id", [f"YB-IC-{i:03d}" for i in range(1, len(supported) + 1)])
    supported.to_csv(run / "screenwide_residual_supported.tsv", sep="\t", index=False)
    summary = dict(n_targets=len(scope), n_match_screenwide=len(full), n_match_targets=full.target_id.nunique(),
                   n_match_compounds=full.inchikey.nunique(), n_both_screenwide=len(supported),
                   n_supported_targets=supported.target_id.nunique(), n_supported_compounds=supported.inchikey.nunique(),
                   n_supported_after_two_variant_bonferroni=int(supported.both_pass_after_two_variant_bonferroni.sum()),
                   n_supported_after_eight_variant_bonferroni=int(supported.both_pass_after_eight_variant_bonferroni.sum()),
                   n_supported_after_twelve_variant_bonferroni=int(supported.both_pass_after_twelve_variant_bonferroni.sum()),
                   variant_count_basis="8 coverage/PC/assay settings plus 4 sparse top-K settings; extra corrections are conservative sensitivity analyses, not biological exclusion rules",
                   match_directional_universe=len(tables["match"]), residual_directional_universe=len(tables["residual"]),
                   directions=supported.direction.value_counts().to_dict(),
                   interpretation="Computational associations. Screen-wide maximum nulls include all pairs, both directions and dose selection in each layer. All match-level associations are retained; this evidence subset is not a biochemical activity validation.")
    (run / "screenwide_export_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return supported


def main():
    run = Path(os.environ.get("YB_SCREENING_RUN", ROOT / "results_a6000/ion_channel_scope_20261012/full_measured_assay"))
    if os.environ.get("YB_SCREENWIDE_EXPORT") == "1":
        export_screenwide_evidence(run)
        return
    if (run / "residual_fullscreen").is_dir():
        selected = export_screenwide_evidence(run)
        out = selected.rename(columns={"candidate_id": "候选编号", "target_id": "靶点",
                                      "inchikey": "化合物_InChIKey", "pubchem_cid": "化合物_CID",
                                      "smiles": "SMILES", "dose": "剂量", "dose_residual": "残差最优剂量"}).copy()
        out.insert(1, "赛道", "人类离子通道候选筛选")
        out["方向"] = out.direction.map({"pos": "条形码相对丰度正向", "neg": "条形码相对丰度负向"})
        out["证据等级"] = "匹配与谱簇残差均通过全库最大统计量检验"
        out["模型与版本"] = MODEL + "；冻结模型完整测量覆盖5865株，HIP/HOP内置换"
        out["备注"] = "计算关联；全部匹配层关联另见screenwide_associations.tsv；非直接人源活性验证"
        columns = ["候选编号", "赛道", "靶点", "化合物_InChIKey", "化合物_CID", "SMILES", "方向",
                   "spearman_rho", "p_match_screenwide", "spearman_rho_residual", "p_residual_screenwide",
                   "剂量", "残差最优剂量", "source_condition_index", "source_condition_index_residual",
                   "证据等级", "模型与版本", "备注",
                   "both_pass_after_two_variant_bonferroni", "both_pass_after_eight_variant_bonferroni",
                   "both_pass_after_twelve_variant_bonferroni"]
        output = (run if os.environ.get("YB_SCREENING_RUN") else ROOT) / "results.csv"
        out[columns].to_csv(output, index=False, encoding="utf-8")
        print(f"{output}：{len(out)}对离子通道计算关联；全库校正P值保持明确字段，不冒充候选内BH q值")
        return
    def result_path(name, legacy):
        return (run / name if run else EH / legacy) / "exec_matrix.tsv"
    pos = pd.read_csv(result_path("confirm_pos", "results_esm2jointdc_confirm100k_pos_20260912"), sep="\t")
    neg = pd.read_csv(result_path("confirm_neg", "results_esm2jointdc_neg_confirm100k_neg_20260912"), sep="\t")
    pos["direction"], neg["direction"] = "+（条形码相对丰度正向）", "-（条形码相对丰度负向）"
    scope_path = run / "target_scope.tsv"
    channel_scope = scope_path.exists()
    if channel_scope:
        import sys
        sys.path.insert(0, str(ROOT / "scripts"))
        from product_execute_hiphop import bh_fdr
        scope = pd.read_csv(scope_path, sep="\t")
        assert set(scope.target_family) == {"ion_channel"}
        assert set(pd.concat([pos, neg]).target_id) <= set(scope.target_id)
        both = pd.concat([pos, neg], ignore_index=True)
        both["q"] = bh_fdr(both.emp_p.to_numpy())
        pos, neg = both.iloc[:len(pos)].copy(), both.iloc[len(pos):].copy()

    # 靶点级证据：残差检验过阈（p<6.25e-4，在册口径）
    resid_sig = set()
    residual_count = sum(len(pd.read_csv(run / "family_panel" / f"resid_allowlist_{d}.tsv", sep="\t"))
                         for d in ("pos", "neg")) if channel_scope else 80
    residual_threshold = .05 / residual_count if residual_count else 0
    for name, d, direction in (("residual_pos", "results_esm2jointdc_resid_pos_20260913", "+（条形码相对丰度正向）"),
                               ("residual_neg", "results_esm2jointdc_resid_neg_neg_20260913", "-（条形码相对丰度负向）")):
        path = result_path(name, d)
        if run and not path.exists():
            continue  # A direction may contain no non-singleton nominated targets.
        r = pd.read_csv(path, sep="\t")
        resid_sig |= {(t, c, direction) for t, c in r[r.emp_p <= residual_threshold][["target_id", "inchikey"]].values}
    # 构造性靶点级：0.9 阈值下无近同谱的单例靶点（在册 build_record）
    build_path = (run / "family_panel" if run else ROOT / "panels/family_specificity_panel_20260913") / "build_record.json"
    build = json.loads(build_path.read_text())
    singles = set(build.get("nominated_singletons", []))

    comp = pd.read_csv(os.environ.get("YB_COMPOUND_TABLE", ROOT / "data" / "compounds_smiles.tsv"), sep="\t")
    comp = comp.loc[comp.inchikey.isin(pd.concat([pos, neg]).inchikey),
                    ["inchikey", "pubchem_cid", "smiles"]].drop_duplicates()
    if comp.inchikey.duplicated().any():
        raise ValueError("conflicting identifiers or structures for a nominated compound")
    cmap = dict(zip(comp.inchikey, comp.smiles))
    cinfo = comp.set_index("inchikey")[["pubchem_cid"]].to_dict("index")

    rows = pd.concat([pos, neg], ignore_index=True).sort_values(["direction", "target_id"])
    out = []
    for i, r in enumerate(rows.itertuples(), 1):
        if (r.target_id, r.inchikey, r.direction) in resid_sig:
            tier = "候选内残差检验过阈"
            note = "去掉家族公共分量后保留关联；候选集合内检验"
        elif r.target_id in singles:
            tier = "任务谱聚类单例"
            note = "0.9 阈值下为单例；未实施家族残差检验"
        else:
            tier = "残差未过阈"
            note = "家族均值关联另见对应检验结果"
        out.append({
            "候选编号": f"YB-{i:03d}", "赛道": "人类离子通道候选筛选" if channel_scope else TRACK, "靶点": r.target_id,
            "化合物_InChIKey": r.inchikey,
            "化合物_CID": cinfo.get(r.inchikey, {}).get("pubchem_cid", ""),
            "SMILES": cmap.get(r.inchikey, ""), "方向": r.direction,
            "spearman_rho": round(r.spearman_rho, 4),
            "emp_p": round(r.emp_p, 6), "q_bh": round(r.q, 8),
            "剂量": r.dose, "证据等级": tier, "模型与版本": MODEL,
            "备注": note,
        })
    df = pd.DataFrame(out)
    output = (run if os.environ.get("YB_SCREENING_RUN") else ROOT) / "results.csv"
    df.to_csv(output, index=False, encoding="utf-8")
    print(f"{output} 生成：{len(df)} 条计算候选；证据分类 {df.证据等级.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
