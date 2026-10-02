#!/usr/bin/env python3
"""Portable entry point for published-result export, a synthetic demo, and full screening."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from scripts.project_io import ROOT, revision, sha256, write_json


def pipeline(mode, output):
    from scripts.check_assets import inventory

    if mode == "full":
        problems = [x for x in inventory() if x["status"] != "ok"]
        if problems:
            raise ValueError(
                "Full run blocked by missing/changed assets:\n"
                + "\n".join(f"  {x['status']}: {x['path']}" for x in problems)
                + "\nSee data/ASSETS.md. No stored confirmation results will be substituted."
            )
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output must be a new or empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    logs.mkdir(exist_ok=True)
    started = time.time()
    record = {
        "mode": mode,
        "code_revision": revision(),
        "status": "running",
        "seed": 42,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy"]
        },
        "steps": [],
    }

    def run(label, script, *args):
        command = [sys.executable, str(ROOT / script), *map(str, args)]
        log = logs / f"{label}.log"
        print(f"[{label}] {log}", flush=True)
        t = time.time()
        with log.open("w") as f:
            proc = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        record["steps"].append(
            {
                "name": label,
                "command": command,
                "returncode": proc.returncode,
                "seconds": round(time.time() - t, 3),
                "log_sha256": sha256(log),
            }
        )
        write_json(output / "run.json", record)
        if proc.returncode:
            raise RuntimeError(f"{label} failed ({proc.returncode}); see {log}")

    try:
        if mode == "demo":
            from scripts.make_demo import design

            t = time.time()
            tasks, response, compounds = design(output / "inputs")
            record["steps"].append(
                {
                    "name": "synthetic-design-and-ridge-fit",
                    "seed": 42,
                    "seconds": round(time.time() - t, 3),
                    "returncode": 0,
                }
            )
            screen_n, confirm_n, family_n, top_n = 99, 199, 199, 3
        else:
            tasks = output / "tasks_ej"
            run(
                "export",
                "scripts/esm2_joint_task_export.py",
                "--config",
                ROOT / "configs/esm2_joint_tasks.json",
                "--output",
                tasks,
            )
            response = ROOT / "raw/tier1_response/strain_response.npz"
            compounds = ROOT / "raw/tier1_response/compounds.tsv.gz"
            screen_n, confirm_n, family_n, top_n = 1000, 100000, 10000, 40
        record["inputs"] = {
            "response": {"path": str(response), "sha256": sha256(response)},
            "compounds": {"path": str(compounds), "sha256": sha256(compounds)},
        }
        dc = output / "tasks_dc"
        run(
            "double_center",
            "scripts/task_axis_double_center.py",
            "--input",
            tasks,
            "--output",
            dc,
        )

        def execute(label, task_dir, sign, n_perm, allowlist=None, rule=None):
            import pandas as pd

            dest = output / label
            if allowlist is not None and pd.read_csv(allowlist, sep="\t").empty:
                dest.mkdir()
                pd.DataFrame(
                    columns=[
                        "target_id",
                        "inchikey",
                        "emp_p",
                        "q",
                        "spearman_rho",
                        "dose",
                    ]
                ).to_csv(dest / "exec_matrix.tsv", sep="\t", index=False)
                write_json(
                    dest / "execute_summary.json",
                    {"status": "skipped-empty-allowlist", "n_pairs": 0},
                )
                record["steps"].append(
                    {
                        "name": label,
                        "returncode": 0,
                        "status": "skipped-empty-allowlist",
                    }
                )
                return
            cfg = {
                "endpoint": "v2-spearman",
                "response_npz": str(response),
                "compound_table": str(compounds),
                "task_dir": str(task_dir),
                "results_dir": str(dest),
                "task_top_k": 6734,
                "n_perm": n_perm,
                "seed": 42,
                "fdr_alpha": 0.1,
                "top_report_per_target": 15,
                "sensitivity_sign": sign,
                "fdr_method": "BH",
                "target_batch_size": 64,
            }
            config = output / "configs" / f"{label}.json"
            write_json(config, cfg)
            args = ["--config", config, "--results-suffix", ""]
            if allowlist is not None:
                args += [
                    "--pair-allowlist",
                    allowlist,
                    "--pair-allowlist-sha256",
                    sha256(allowlist),
                    "--allowlist-two-stage-rule-doc",
                    rule,
                    "--allowlist-selection-basis",
                    "Outcome-derived two-stage nomination; exploratory family-only FDR; see supplied rule document.",
                ]
            run(label, "scripts/product_execute_hiphop.py", *args)

        for name, sign in [("pos", 1), ("neg", -1)]:
            execute("screen_" + name, dc, sign, screen_n)
        panel = output / "confirmation_panel"
        run(
            "nominate",
            "scripts/ej_dc_nominate.py",
            "--pos-screen",
            output / "screen_pos",
            "--neg-screen",
            output / "screen_neg",
            "--top-n",
            top_n,
            "--screen-permutations",
            screen_n,
            "--out",
            panel,
        )
        rule = ROOT / (
            "data/demo/README.md"
            if mode == "demo"
            else "docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md"
        )
        for name, sign in [("pos", 1), ("neg", -1)]:
            execute(
                "confirm_" + name,
                dc,
                sign,
                confirm_n,
                panel / f"pair_allowlist_{name}40.tsv",
                rule,
            )
        family = output / "family"
        run(
            "family_build",
            "scripts/family_specificity_build.py",
            "--tasks",
            dc,
            "--confirmation-panel",
            panel,
            "--output",
            family,
        )
        for kind, folder in [("fam", "tasks_fammean"), ("resid", "tasks_resid")]:
            for name, sign in [("pos", 1), ("neg", -1)]:
                execute(
                    f"{kind}_{name}",
                    family / folder,
                    sign,
                    family_n,
                    family / "panel" / f"{kind}_allowlist_{name}.tsv",
                    ROOT
                    / (
                        "data/demo/README.md"
                        if mode == "demo"
                        else "docs/cross_species_match/FAMILY_SPECIFICITY_PROTOCOL.md"
                    ),
                )
        args = [
            "--run-dir",
            output,
            "--compounds",
            compounds,
            "--output",
            output / "results.csv",
        ]
        if mode == "demo":
            args += ["--demo"]
        run("predict", "predict.py", *args)
        record["status"] = "complete"
        record["result_sha256"] = sha256(output / "results.csv")
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = str(exc)
        raise
    finally:
        record["seconds"] = round(time.time() - started, 3)
        write_json(output / "run.json", record)
    print(f'Completed {mode}: {output / "results.csv"}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mode", choices=["export", "demo", "full", "check"], default="export"
    )
    p.add_argument(
        "--output", type=Path, help="CSV file for export; new directory for demo/full"
    )
    args = p.parse_args()
    if args.mode == "export":
        from predict import export

        export(output=args.output)
    elif args.mode == "check":
        from scripts.check_assets import inventory

        items = inventory()
        print(json.dumps(items, ensure_ascii=False, indent=2))
        return 0 if all(x["status"] == "ok" for x in items) else 2
    else:
        pipeline(args.mode, (args.output or ROOT / "outputs" / args.mode).resolve())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2)
