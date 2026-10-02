#!/usr/bin/env python3
"""Build a small review bundle, or a full bundle only after asset/training checks pass."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, revision
from scripts.check_assets import inventory


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile", choices=["review", "full"], default="review")
    p.add_argument("--output", type=Path, default=ROOT / "dist/YeastBridge-review.zip")
    args = p.parse_args()
    tracked = (
        subprocess.check_output(
            ["git", "ls-files", "--cached", "-z"],
            cwd=ROOT,
        )
        .decode()
        .split("\0")
    )
    paths = {x for x in tracked if x and (ROOT / x).is_file()}
    paths = {
        x
        for x in paths
        if not x.startswith(("outputs/", "dist/", ".venv/")) and "__pycache__" not in x
    }
    if args.profile == "full":
        issues = [x for x in inventory() if x["status"] != "ok"]
        if issues:
            p.error(
                "Full bundle blocked: missing or changed assets. Run main.py --mode check."
            )
        manifest = ROOT / "models/training_manifest.json"
        if not manifest.is_file():
            p.error(
                "Full bundle blocked: missing B2 training_manifest.json; see models/TRAINING.md."
            )
        evidence = json.loads(manifest.read_text())
        paths.add("models/training_manifest.json")
        for key in ["entrypoint", "environment", "logs", "splits"]:
            values = evidence.get(key, [])
            if isinstance(values, str):
                values = [values]
            if not values or any(
                not (ROOT / v).is_file()
                or not (ROOT / v).resolve().is_relative_to(ROOT)
                for v in values
            ):
                p.error(f"Full bundle blocked: missing B2 training evidence: {key}")
            paths.update(values)
        paths.update(
            [
                "raw/tier1_models/route_b/final_model.pt",
                "raw/tier1_response/strain_response.npz",
            ]
        )
    else:
        roots = {
            "README.md",
            "CITATION.cff",
            "main.py",
            "predict.py",
            "run.sh",
            "requirements.txt",
            "requirements-full.txt",
            "requirements.lock.txt",
            "requirements-full.lock.txt",
            "NOTICE.md",
            "MANIFEST.sha256",
            ".gitattributes",
            ".gitignore",
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
        )
        keep_results = {
            "results.csv",
            *[
                f"results_a6000/execute_hiphop/{name}/exec_matrix.tsv"
                for name in [
                    "results_esm2jointdc_confirm100k_pos_20260912",
                    "results_esm2jointdc_neg_confirm100k_neg_20260912",
                    "results_esm2jointdc_resid_pos_20260913",
                    "results_esm2jointdc_resid_neg_neg_20260913",
                ]
            ],
        }
        paths = {
            x
            for x in paths
            if x in roots or x in keep_results or x.startswith(prefixes)
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sums = []
    with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name in sorted(paths):
            data = (ROOT / name).read_bytes()
            z.writestr("YeastBridge/" + name, data)
            sums.append(hashlib.sha256(data).hexdigest() + "  " + name)
        z.writestr("YeastBridge/MANIFEST.release.sha256", "\n".join(sums) + "\n")
        z.writestr(
            "YeastBridge/BUNDLE.json",
            json.dumps(
                {
                    "profile": args.profile,
                    "code_revision": revision(),
                    "files": len(paths),
                    "purpose": (
                        "published export and synthetic demo; not a complete research submission"
                        if args.profile == "review"
                        else "full assets and supplied training evidence; scientific/licensing review remains required"
                    ),
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
