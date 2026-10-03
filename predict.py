#!/usr/bin/env python3
"""Export the published candidate table, or results from an explicitly selected run."""
import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.project_io import ROOT, revision, sha256, write_json

TRACK = "赛道三（离子通道、GPCR 等药靶的小分子药物虚拟筛选）"
MODEL = "YeastBridge EJ-dc 任务轴（RRF20: direct_pc1 + aligned_C + b2_verbatim，共模双去除）"
HISTORICAL = {
    "confirm_pos": "results_esm2jointdc_confirm100k_pos_20260912",
    "confirm_neg": "results_esm2jointdc_neg_confirm100k_neg_20260912",
    "resid_pos": "results_esm2jointdc_resid_pos_20260913",
    "resid_neg": "results_esm2jointdc_resid_neg_neg_20260913",
}


def export(run_dir=None, output=None, compounds=None, demo=False, legacy=False):
    published = run_dir is None
    submission_path = ROOT / "models/submission.json"
    if published and not legacy and submission_path.is_file():
        submission = json.loads(submission_path.read_text())
        for name, expected in submission["export_sources"].items():
            if sha256(ROOT / name) != expected:
                raise ValueError(f"Published source checksum mismatch: {name}")
        run_dir = ROOT / submission["run_dir"]
    historical = run_dir is None
    base = (
        ROOT / "results_a6000/execute_hiphop" if historical else Path(run_dir).resolve()
    )
    output = Path(output or ROOT / "results.csv").resolve()
    compound_path = Path(compounds or ROOT / "data/compounds_smiles.tsv").resolve()
    names = HISTORICAL if historical else {k: k for k in HISTORICAL}
    sources = []

    def read(name):
        p = base / names[name] / "exec_matrix.tsv"
        sources.append(p)
        frame = pd.read_csv(p, sep="\t")
        required = {"target_id", "inchikey", "emp_p"}
        if not required.issubset(frame):
            raise ValueError(f"{p}: missing columns {sorted(required - set(frame))}")
        return frame

    pos, neg = read("confirm_pos"), read("confirm_neg")
    pos["direction"], neg["direction"] = "+z（耐药富集）", "-z（超敏）"
    # Keep published statistical values. New runs use BH over BOTH directions.
    rows = pd.concat([pos, neg], ignore_index=True)
    if rows.empty:
        raise ValueError("No confirmed candidate records to export")
    if not historical:
        from scripts.product_execute_hiphop import bh_fdr

        rows["q"] = bh_fdr(rows.emp_p.to_numpy())
    threshold = 6.25e-4 if historical else 0.05 / len(rows)
    resid_sig = set()
    for name in ("resid_pos", "resid_neg"):
        frame = read(name)
        resid_sig |= set(
            map(
                tuple,
                frame.loc[frame.emp_p < threshold, ["target_id", "inchikey"]].values,
            )
        )
    build_path = (
        ROOT / "panels/family_specificity_panel_20260913/build_record.json"
        if historical
        else base / "family/panel/build_record.json"
    )
    build = json.loads(build_path.read_text())
    singles = set(build["nominated_singletons"])
    comp = pd.read_csv(compound_path, sep="\t").fillna("")
    # The source table contains repeated records from multiple screens. Collapse only
    # consistent metadata; never arbitrarily choose among conflicting structures.
    cmap = {}
    for key in set(rows.inchikey):
        group = comp[comp.inchikey == key]
        if group.empty:
            continue
        entry = {}
        for column in ("smiles", "pubchem_cid"):
            values = (
                list(dict.fromkeys(v for v in group[column] if v != ""))
                if column in group
                else []
            )
            if len(values) > 1:
                raise ValueError(
                    f"Conflicting {column} metadata for {key}; resolve before export"
                )
            entry[column] = values[0] if values else ""
        cmap[key] = entry
    missing = set(rows.inchikey) - set(cmap)
    if missing:
        raise ValueError(f"Missing compound metadata: {sorted(missing)}")
    if any(not cmap[k].get("smiles") for k in set(rows.inchikey)):
        raise ValueError(
            "A nominated compound has no SMILES; resolve metadata before export"
        )
    unit_map = {}
    if not demo:
        condition_path = ROOT / "data/response_conditions.tsv"
        conditions = pd.read_csv(condition_path, sep="\t")
        sources.append(condition_path)
        for condition in conditions.itertuples():
            key = (condition.inchikey, float(condition.dose))
            unit_map.setdefault(key, set()).add(condition.dose_unit)
    rows = rows.sort_values(["direction", "target_id"], kind="stable")
    model_version = "EJ-dc-B2-20260912"
    run_record = base / "run.json"
    if not historical and not demo and run_record.is_file():
        model_version = (
            json.loads(run_record.read_text())
            .get("model", {})
            .get("version", model_version)
        )
        sources.append(run_record)
    out = []
    for i, r in enumerate(rows.itertuples(), 1):
        if not historical and r.q >= 0.1:
            tier, note = "候选（未通过确认阈值）", "联合确认家族 BH q≥0.1"
        elif (r.target_id, r.inchikey) in resid_sig:
            tier, note = "靶点级（残差检验过阈）", "去掉家族公共分量后仍显著"
        elif r.target_id in singles:
            tier, note = "靶点级（构造性，无近同谱）", "0.9 阈值下为单例靶点"
        else:
            tier = "家族级" if historical else "确认关联（家族内未分辨）"
            note = (
                "家族均值轴复现（12/12 在册）"
                if historical
                else "残差未通过升级阈值；家族检验见本次运行日志"
            )
        if demo:
            tier, note = "合成演示，非科研候选", "合成表示与响应，仅验证计算流程"
        units = (
            {"arbitrary demo units"}
            if demo
            else unit_map.get((r.inchikey, float(r.dose)), set())
        )
        if len(units) != 1:
            raise ValueError(
                f"Missing or ambiguous dose unit: {r.inchikey}, {r.dose}: {units}"
            )
        out.append(
            {
                "候选编号": f'{"DEMO" if demo else "YB"}-{i:03d}',
                "赛道": TRACK,
                "靶点": r.target_id,
                "化合物_InChIKey": r.inchikey,
                "化合物_CID": cmap[r.inchikey].get("pubchem_cid", ""),
                "SMILES": cmap[r.inchikey]["smiles"],
                "方向": r.direction,
                "spearman_rho": round(r.spearman_rho, 4),
                "emp_p": round(r.emp_p, 6),
                "q_bh": round(r.q, 8),
                "剂量": r.dose,
                "剂量单位": next(iter(units)),
                "证据等级": tier,
                "模型与版本": (
                    "synthetic-demo-v1 (not ESM-2/B2 weights)"
                    if demo
                    else (MODEL if historical else MODEL + "; " + model_version)
                ),
                "备注": note,
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(out).to_csv(output, index=False, encoding="utf-8")
    sources += [build_path, compound_path]
    write_json(
        output.with_suffix(".metadata.json"),
        {
            "mode": (
                "published-export"
                if published
                else "synthetic-demo" if demo else "new-run"
            ),
            "code_revision": revision(),
            "rows": len(out),
            "output_sha256": sha256(output),
            "sorting": "direction then target symbol; candidate IDs are not a potency ranking",
            "dose_unit": (
                "arbitrary demo units"
                if demo
                else "per-row 剂量单位 column, recovered from the frozen response file"
            ),
            "q_scope": (
                "per-direction confirmation family (historical values retained)"
                if historical
                else "union of both directional confirmation families"
            ),
            "selection_scope": "BH does not account for stage-one nomination; exploratory associations",
            "sources": {
                str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha256(
                    p
                )
                for p in sources
            },
        },
    )
    print(f"{output}: {len(out)} candidate records")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Read a new main.py run instead of published records",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "results.csv")
    parser.add_argument("--compounds", type=Path)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument(
        "--legacy", action="store_true", help="Export archived September records"
    )
    args = parser.parse_args()
    if args.demo and not args.run_dir:
        parser.error("--demo requires --run-dir")
    export(args.run_dir, args.output, args.compounds, args.demo, args.legacy)


if __name__ == "__main__":
    main()
