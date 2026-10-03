#!/usr/bin/env python3
"""Validate that the submitted model, training evidence and candidate table agree."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256, write_json
from scripts.check_assets import inventory


def read(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError(f"Path outside repository: {relative}")
    return json.loads(path.read_text())


def verify_training():
    problems = []
    try:
        manifest = read("models/training_manifest.json")
        if manifest.get("final_checkpoint_epoch_logs_verified") is not True:
            return ["matching final-checkpoint epoch logs have not been supplied"]
        run = read(manifest["run_record"])
        config = read(manifest["config"])
        metrics = read(manifest["epoch_metrics"])
        meta = read(manifest["train_meta"])
        if (
            run["status"] != "complete"
            or run.get("returncode") != 0
            or run["mode"] != "training"
        ):
            problems.append("training did not complete successfully")
        if len(metrics) != config["epochs"] or [m["epoch"] for m in metrics] != list(
            range(1, config["epochs"] + 1)
        ):
            problems.append("missing or out-of-order epoch metrics")
        for row in metrics:
            for key in ["train_mse_mask_weighted", "validation_mse_mask_weighted"]:
                if not math.isfinite(row[key]) or row[key] < 0:
                    problems.append(f"invalid metric: epoch {row['epoch']} {key}")
        checkpoint_sha = sha256(ROOT / manifest["final_checkpoint_path"])
        if (
            not checkpoint_sha
            == manifest["final_checkpoint_sha256"]
            == meta["checkpoint_sha256"]
            == run["artifacts"]["model/final_model.pt"]
        ):
            problems.append("checkpoint does not match training outputs")
        for relative, artifact in manifest["evidence_artifacts"].items():
            if sha256(ROOT / relative) != run["artifacts"][artifact]:
                problems.append(f"training evidence checksum mismatch: {relative}")
        for key, path in [
            ("source_sha256", "scripts/scf_finetune_yeast.py"),
            ("wrapper_sha256", "train.py"),
        ]:
            if sha256(ROOT / path) != run[key]:
                problems.append(f"training source differs from recorded run: {path}")
        log = (ROOT / manifest["logs"][0]).read_text()
        if any(
            f"[ep{epoch}] train=" not in log for epoch in range(1, config["epochs"] + 1)
        ):
            problems.append("stdout lacks complete epoch train/validation records")
        if sha256(ROOT / manifest["config"]) != run["config_sha256"]:
            problems.append("training configuration differs from recorded run")
        if run["required_assets"] != config["required_assets"]:
            problems.append("training input versions disagree")
        if (
            meta["implementation"] != config["implementation"]
            or run["implementation"] != config["implementation"]
        ):
            problems.append("training implementation versions disagree")
    except (OSError, KeyError, ValueError, TypeError) as error:
        problems.append(f"training evidence unavailable or invalid: {error}")
    return problems


def verify_submission():
    checks = {}
    checks["training"] = verify_training()
    try:
        checks["runtime_assets"] = [
            f"{x['path']}: {x['status']}" for x in inventory() if x["status"] != "ok"
        ]
        submission = read("models/submission.json")
        run = read(submission["run_dir"] + "/run.json")
        manifest = read("models/training_manifest.json")
        training = read(manifest["run_record"])
        task = read(submission["task_config"])
        assets = read(task["method"]["v7_config"])
        checks["run"] = []
        if (
            run["status"] != "complete"
            or run["mode"] != "full"
            or any(s.get("returncode") != 0 for s in run["steps"])
        ):
            checks["run"].append("full screening did not complete successfully")
        expected_steps = [
            "export",
            "double_center",
            "screen_pos",
            "screen_neg",
            "nominate",
            "confirm_pos",
            "confirm_neg",
            "family_build",
            "fam_pos",
            "fam_neg",
            "resid_pos",
            "resid_neg",
            "uncertainty",
            "predict",
        ]
        if [step["name"] for step in run["steps"]] != expected_steps:
            checks["run"].append("full screening step sequence is incomplete")
        if run["model"]["checkpoint_sha256"] != manifest["final_checkpoint_sha256"]:
            checks["run"].append("candidate pipeline used a different checkpoint")
        if not (
            run["model"]["version"]
            == submission["model_version"]
            == task["model_version"]
            and run["model"]["task_config_sha256"]
            == sha256(ROOT / submission["task_config"])
            and run["model"]["checkpoint_sha256"]
            == assets["inputs"]["route_b_model"]["sha256"]
            and run["model"]["gene_table_sha256"]
            == assets["inputs"]["route_b_table"]["sha256"]
            == training["artifacts"]["model/gene_table_final.npy"]
        ):
            checks["run"].append(
                "published configuration differs from the completed run"
            )
        if run["result_sha256"] != sha256(ROOT / "results.csv"):
            checks["run"].append("published candidates differ from the full-run result")
        for step in run["steps"]:
            if "log_sha256" in step:
                log = ROOT / submission["run_dir"] / "logs" / (step["name"] + ".log")
                if sha256(log) != step["log_sha256"]:
                    checks["run"].append(f"run log checksum mismatch: {step['name']}")
        from predict import export

        with tempfile.TemporaryDirectory() as temp:
            out = export(output=Path(temp) / "results.csv")
            checks["export"] = (
                []
                if out.read_bytes() == (ROOT / "results.csv").read_bytes()
                else ["published results differ from regenerated export"]
            )
        with (ROOT / "results.csv").open() as file:
            rows = list(csv.DictReader(file))
        required = [
            "候选编号",
            "赛道",
            "SMILES",
            "靶点",
            "化合物_InChIKey",
            "模型与版本",
            "备注",
            "剂量单位",
        ]
        checks["schema"] = []
        if len(rows) != submission["rows"] or len({r["候选编号"] for r in rows}) != len(
            rows
        ):
            checks["schema"].append("candidate count or unique IDs disagree")
        for row in rows:
            if any(not row.get(key) for key in required):
                checks["schema"].append(
                    f"required field missing: {row.get('候选编号')}"
                )
            if submission["model_version"] not in row["模型与版本"]:
                checks["schema"].append(
                    "result lacks the exact submitted model version"
                )
            if (
                not -1 <= float(row["spearman_rho"]) <= 1
                or not 0 < float(row["emp_p"]) <= 1
                or not 0 <= float(row["q_bh"]) <= 1
            ):
                checks["schema"].append("invalid statistical value")
        datasets = read("data/datasets.json")
        checks["data_disclosure"] = []
        for item in datasets["datasets"]:
            if any(
                not item.get(key)
                for key in [
                    "name",
                    "source",
                    "content_version",
                    "role",
                    "acquired_for_submission",
                    "terms_url",
                ]
            ):
                checks["data_disclosure"].append(
                    f"missing disclosure: {item.get('name')}"
                )
    except (OSError, KeyError, ValueError, TypeError) as error:
        checks["submission"] = [str(error)]
    report = {
        "status": "pass" if all(not issues for issues in checks.values()) else "fail",
        "checks": checks,
        "scope": "Attachment 5 artifact consistency and default CSV contract; not an organizer acceptance decision",
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = verify_submission()
    if args.output:
        write_json(args.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
