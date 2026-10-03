#!/usr/bin/env python3
"""Create a review, inference or complete reproducibility submission bundle."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, revision, sha256
from scripts.check_assets import inventory


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--profile", choices=["review", "inference", "full"], default="review"
    )
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    args.output = args.output or ROOT / f"dist/YeastBridge-{args.profile}.zip"
    tracked = (
        subprocess.check_output(["git", "ls-files", "--cached", "-z"], cwd=ROOT)
        .decode()
        .split("\0")
    )
    roots = {
        "README.md",
        "CITATION.cff",
        "main.py",
        "predict.py",
        "train.py",
        "run.sh",
        "requirements-training.txt",
        "requirements-embeddings.txt",
        "requirements-training.lock.txt",
        "environment-preprocessing.yml",
        "requirements.txt",
        "requirements-full.txt",
        "requirements.lock.txt",
        "requirements-full.lock.txt",
        "NOTICE.md",
        "MANIFEST.sha256",
        ".gitattributes",
        ".gitignore",
        "results.csv",
    }
    prefixes = (
        "scripts/",
        "tests/",
        "data/",
        "models/",
        "docs/",
        "configs/",
        "configs_a6000_frozen/",
        "panels/",
        "repro/",
    )
    historical = {
        f"results_a6000/execute_hiphop/{name}/exec_matrix.tsv"
        for name in [
            "results_esm2jointdc_confirm100k_pos_20260912",
            "results_esm2jointdc_neg_confirm100k_neg_20260912",
            "results_esm2jointdc_resid_pos_20260913",
            "results_esm2jointdc_resid_neg_neg_20260913",
        ]
    }
    model_call_records = {
        "raw/tier1_esm2/human_810/run_info.json",
        "raw/tier1_esm2/human_810/build_manifest.tsv",
        "raw/tier1_esm2/yeast_650m/run_info.json",
    }
    paths = {
        name
        for name in tracked
        if name
        and (ROOT / name).is_file()
        and (
            name in roots
            or name in historical
            or name in model_call_records
            or name.startswith(prefixes)
        )
        and "__pycache__" not in name
    }
    if args.profile in {"inference", "full"}:
        items = inventory()
        if any(item["status"] != "ok" for item in items):
            p.error(
                "Bundle blocked: missing or changed runtime assets; run main.py --mode check"
            )
        paths.update(item["path"] for item in items)
    if args.profile == "full":
        from scripts.verify_submission import verify_submission

        report = verify_submission()
        issues = [
            f"{name}: {issue}"
            for name, group in report["checks"].items()
            for issue in group
        ]
        if issues:
            p.error("Full bundle blocked: " + "; ".join(issues))
        archive = "raw/training/b2-inputs.tar.gz"
        asset = json.loads((ROOT / "data/asset_downloads.json").read_text())["assets"][
            "training_inputs"
        ]
        if not (ROOT / archive).is_file() or sha256(ROOT / archive) != asset["sha256"]:
            p.error(
                "Full bundle requires the verified training-input archive; run install_training_assets.py"
            )
        paths.add(archive)
        paths.add("raw/externalvalidation/mappings/gene_master.tsv")
    # Every package uses one explicit allowlist; historical task matrices and
    # evaluation-only model assets are never pulled in by a broad raw/ prefix.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sums = []
    with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name in sorted(paths):
            z.write(ROOT / name, "YeastBridge/" + name)
            sums.append(sha256(ROOT / name) + "  " + name)
        z.writestr("YeastBridge/MANIFEST.release.sha256", "\n".join(sums) + "\n")
        z.writestr(
            "YeastBridge/BUNDLE.json",
            json.dumps(
                {
                    "profile": args.profile,
                    "code_revision": revision(),
                    "files": len(paths),
                    "purpose": {
                        "review": "Source, documentation, published export and synthetic demo",
                        "inference": "Source and all final-pipeline runtime assets",
                        "full": "Runtime assets, final checkpoint, prepared training inputs and matching training evidence; obtain the upstream initialization checkpoint from its official source",
                    }[args.profile],
                },
                indent=2,
            )
            + "\n",
        )
    print(
        f"{args.output}: {len(paths)} files, {args.output.stat().st_size/1048576:.2f} MiB ({args.profile})"
    )


if __name__ == "__main__":
    main()
