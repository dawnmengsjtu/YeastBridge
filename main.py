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
from datetime import datetime, timezone

from scripts.project_io import ROOT, revision, sha256, write_json


def pipeline(mode, output, threads=4, task_config=None):
    from scripts.check_assets import inventory, default_task_config

    # Apply before importing NumPy and inherit the same limits in every worker.
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = str(threads)
    if mode == "full":
        task_config = Path(task_config or default_task_config()).resolve()
        problems = [
            x for x in inventory(task_config=task_config) if x["status"] != "ok"
        ]
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
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": {
            str(p.relative_to(ROOT)): sha256(p)
            for p in [
                ROOT / "main.py",
                ROOT / "predict.py",
                *sorted((ROOT / "scripts").glob("*.py")),
            ]
        },
        "code_revision": revision(),
        "status": "running",
        "seed": 42,
        "blas_threads": threads,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            p: importlib.metadata.version(p)
            for p in (
                ["numpy", "pandas", "scipy", "torch", "anndata", "h5py", "fair-esm"]
                if mode == "full"
                else ["numpy", "pandas", "scipy"]
            )
        },
        "steps": [],
    }

    if mode == "full":
        task = json.loads(task_config.read_text())
        v7 = json.loads((ROOT / task["method"]["v7_config"]).read_text())
        record["model"] = {
            "version": task.get("model_version", "EJ-dc-B2-20260912"),
            "task_config_sha256": sha256(task_config),
            "checkpoint_sha256": v7["inputs"]["route_b_model"]["sha256"],
            "gene_table_sha256": v7["inputs"]["route_b_table"]["sha256"],
        }
        write_json(output / "configs/task.json", task)
        write_json(output / "configs/model.json", v7)

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
                task_config,
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
        run("uncertainty", "scripts/assess_uncertainty.py", "--run-dir", output)
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
        record["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(output / "run.json", record)
    # predict reads run.json for the model version before this wrapper finalizes it.
    # Bind its metadata to the completed record, rather than the transient hash.
    if mode == "full":
        metadata_path = output / "results.metadata.json"
        metadata = json.loads(metadata_path.read_text())
        run_path = output / "run.json"
        key = (
            str(run_path.relative_to(ROOT))
            if run_path.is_relative_to(ROOT)
            else str(run_path)
        )
        metadata["sources"][key] = sha256(run_path)
        write_json(metadata_path, metadata)
    print(f'Completed {mode}: {output / "results.csv"}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mode", choices=["export", "demo", "full", "check"], default="export"
    )
    p.add_argument(
        "--output", type=Path, help="CSV file for export; new directory for demo/full"
    )
    p.add_argument(
        "--threads",
        type=int,
        default=4,
        help="BLAS/OpenMP thread limit for demo/full (default: 4)",
    )
    p.add_argument(
        "--task-config",
        type=Path,
        help="Explicit task/model configuration for full/check; defaults to submission version",
    )
    p.add_argument(
        "--legacy",
        action="store_true",
        help="Export the archived September candidate table",
    )
    args = p.parse_args()
    if args.threads < 1:
        p.error("--threads must be positive")
    if args.mode == "export":
        from predict import export

        export(output=args.output, legacy=args.legacy)
    elif args.mode == "check":
        from scripts.check_assets import inventory, default_task_config

        items = inventory(task_config=args.task_config)
        print(json.dumps(items, ensure_ascii=False, indent=2))
        return 0 if all(x["status"] == "ok" for x in items) else 2
    else:
        pipeline(
            args.mode,
            (args.output or ROOT / "outputs" / args.mode).resolve(),
            args.threads,
            args.task_config,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2)
