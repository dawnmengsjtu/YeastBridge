#!/usr/bin/env python3
"""Connect a completed B2 training run to the full candidate-screening entry point."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256, write_json

EXPORT_INPUTS = {
    "esm2_human",
    "route_b_table",
    "route_b_model",
    "route_b_gene_order",
    "yeast_esm2",
    "yeast_metadata",
    "v51_pair_list",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--version", help="Human-readable model version; defaults to checkpoint digest"
    )
    args = parser.parse_args()
    training = args.training_run.resolve()
    output = args.output.resolve()
    if not training.is_relative_to(ROOT) or not output.is_relative_to(ROOT):
        parser.error(
            "Keep training and configuration directories inside the repository for portable paths"
        )
    if output.exists() and any(output.iterdir()):
        parser.error("Output must be a new or empty directory")
    run = json.loads((training / "run.json").read_text())
    if (
        run.get("status") != "complete"
        or run.get("mode") != "training"
        or run.get("returncode") != 0
    ):
        parser.error(
            "A successfully completed full training run is required; smoke output is not a model"
        )
    if run.get("implementation") != "submission-v2":
        parser.error("This entry expects the submission-v2 training recipe")
    for name in [
        "model/final_model.pt",
        "model/gene_table_final.npy",
        "model/split.tsv",
        "model/epoch_metrics.json",
    ]:
        if sha256(training / name) != run["artifacts"][name]:
            parser.error(f"Training output checksum mismatch: {name}")
    config = json.loads((training / "config.json").read_text())
    expected = json.loads((ROOT / "configs/b2_training_v2.json").read_text())
    if config["required_assets"] != expected["required_assets"]:
        parser.error(
            "The training input/gene-order version differs from this screening configuration"
        )
    legacy = json.loads((ROOT / "configs/v7_b2_joint_formal.json").read_text())
    inputs = {k: v for k, v in legacy["inputs"].items() if k in EXPORT_INPUTS}
    for key, name in [
        ("route_b_model", "final_model.pt"),
        ("route_b_table", "gene_table_final.npy"),
    ]:
        path = training / "model" / name
        inputs[key] = {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
    pairs = inputs["v51_pair_list"]
    pairs["sha256"] = sha256(ROOT / pairs["path"])
    task = json.loads((ROOT / "configs/esm2_joint_tasks.json").read_text())
    task.pop("user_rulings", None)
    task["model_version"] = (
        args.version or "B2-submission-v2-" + inputs["route_b_model"]["sha256"][:12]
    )
    task["method"]["setup_mode"] = "task-export-only"
    task["method"][
        "train_pairs"
    ] = "Frozen project split: 753 pairs; 510 training pairs with available representations. No competition-provided training split is implied."
    task["method"]["v7_config"] = str((output / "assets.json").relative_to(ROOT))
    task["training_run_sha256"] = sha256(training / "run.json")
    task["note"] = (
        "Original EJ hyperparameters retained; B2 checkpoint replaced by the explicitly recorded new training run. Historical benchmark scores are not reassigned to this model."
    )
    write_json(output / "assets.json", {"stage_root": ".", "inputs": inputs})
    write_json(output / "task.json", task)
    print(f"Model configuration: {(output/'task.json').relative_to(ROOT)}")
    print(
        f"Run: python main.py --mode full --task-config {(output/'task.json').relative_to(ROOT)} --output outputs/b2-screen"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, KeyError, ValueError) as error:
        print(f"Model configuration failed: {error}", file=sys.stderr)
        raise SystemExit(2)
