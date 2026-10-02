#!/usr/bin/env python3
"""Rebuild HIP/HOP responses from frozen expression or the 18 official CEL archives."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256
from scripts.download_asset import download


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--from-cel",
        action="store_true",
        help="Download ~2.8 GB of official CEL archives and run R preprocessing",
    )
    p.add_argument("--rscript", default="Rscript")
    p.add_argument(
        "--expression",
        type=Path,
        help="Use a local expr_full.tsv.gz instead of downloading",
    )
    args = p.parse_args()
    if args.from_cel and args.expression:
        p.error("--from-cel and --expression are mutually exclusive")
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        p.error("Output must be a new or empty directory")
    out.mkdir(parents=True, exist_ok=True)
    if args.from_cel:
        cache = ROOT / "raw/preprocessing"
        cache.mkdir(parents=True, exist_ok=True)
        manifest = json.loads(
            (ROOT / "data/provenance/hiphop/hiphop_readiness_final.json").read_text()
        )
        arrays = {}
        cel_dir = cache / "cel"
        cel_dir.mkdir(exist_ok=True)
        for item in manifest["archives"]["records"]:
            path = cache / Path(item["path"]).name
            if not path.exists() or sha256(path) != item["sha256"]:
                temp = path.with_suffix(".partial")
                print(f"Downloading {path.name}", flush=True)
                with urllib.request.urlopen(item["url"], timeout=120) as src, temp.open(
                    "wb"
                ) as dst:
                    shutil.copyfileobj(src, dst)
                if sha256(temp) != item["sha256"]:
                    temp.unlink()
                    raise ValueError(
                        f"Official CEL archive checksum mismatch: {path.name}"
                    )
                temp.replace(path)
            with zipfile.ZipFile(path) as z:
                for member in z.infolist():
                    if member.is_dir() or not member.filename.upper().endswith(".CEL"):
                        continue
                    name = Path(member.filename).name
                    if name in arrays:
                        raise ValueError(f"Duplicate CEL name in archives: {name}")
                    target = cel_dir / name
                    with z.open(member) as src, target.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                    arrays[name] = target
        # Preserve the frozen column order, rather than filesystem/ZIP iteration order.
        with (ROOT / "data/response_conditions.tsv").open() as f:
            order = [row["array_name"] for row in csv.DictReader(f, delimiter="\t")]
        if set(order) != set(arrays) or len(order) != len(arrays):
            raise ValueError(
                "CEL array inventory differs from the frozen response conditions"
            )
        cel_list = out / "cel_list.tsv"
        with cel_list.open("w") as f:
            for name in order:
                f.write(f"{name}\t{arrays[name]}\n")
        expression = out / "expr_full.tsv.gz"
        command = [
            args.rscript,
            str(ROOT / "scripts/preprocessing/hiphop_cel_process.R"),
            str(cel_list),
            str(ROOT / "data/provenance/hiphop/probeset_map_xy.tsv"),
            str(expression),
            str(out / "cel_summary.txt"),
        ]
        with (out / "cel.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        expression = args.expression or download("hiphop_expression")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/preprocess_hiphop.py"),
            "--expr",
            str(expression.resolve()),
            "--sdrf",
            str(ROOT / "data/provenance/hiphop/E-MTAB-2391.sdrf.txt"),
            "--output-npz",
            str(out / "strain_response.npz"),
            "--output-json",
            str(out / "record.json"),
        ],
        check=True,
    )
    print(f"Rebuilt response and provenance: {out}")


if __name__ == "__main__":
    main()
