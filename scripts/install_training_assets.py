#!/usr/bin/env python3
"""Install frozen B2 inputs and a hash-checked upstream checkpoint for training."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256
from scripts.download_asset import download


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--inputs",
        type=Path,
        help="Local b2-inputs.tar.gz; otherwise download the release asset",
    )
    p.add_argument(
        "--checkpoint",
        type=Path,
        help="Official scFoundation models.ckpt (see models/TRAINING.md)",
    )
    args = p.parse_args()
    archive = args.inputs or download("training_inputs")
    expected = json.loads((ROOT / "data/asset_downloads.json").read_text())["assets"][
        "training_inputs"
    ]["sha256"]
    if sha256(archive) != expected:
        p.error("Training input archive checksum mismatch")
    destination = ROOT / "raw/training"
    destination.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((ROOT / "configs/b2_training.json").read_text())
    expected_files = {
        Path(cfg["required_assets"][k]["path"]).name: cfg["required_assets"][k][
            "sha256"
        ]
        for k in ["counts", "protein_matrix"]
    }
    with tempfile.TemporaryDirectory(dir=destination) as temp:
        with tarfile.open(archive, "r:gz") as source:
            members = source.getmembers()
            if len(members) != len(expected_files) or {m.name for m in members} != set(
                expected_files
            ):
                p.error("Unexpected members in training archive")
            for member in members:
                if not member.isfile():
                    p.error("Training archive must contain regular files only")
                target = Path(temp) / member.name
                with source.extractfile(member) as src, target.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
                if sha256(target) != expected_files[member.name]:
                    p.error(f"Training member checksum mismatch: {member.name}")
        for name in expected_files:
            shutil.move(str(Path(temp) / name), destination / name)
    shutil.copyfile(
        ROOT / "models/training/evidence/corpus_cells.tsv",
        destination / "corpus_cells.tsv",
    )
    if args.checkpoint:
        expected_checkpoint = cfg["required_assets"]["pretrained_checkpoint"]["sha256"]
        if sha256(args.checkpoint) != expected_checkpoint:
            p.error("Pretrained checkpoint checksum mismatch; no model was installed")
        target = ROOT / cfg["required_assets"]["pretrained_checkpoint"]["path"]
        if args.checkpoint.resolve() != target.resolve():
            shutil.copyfile(args.checkpoint, target)
    print(
        "Training inputs verified. The upstream models.ckpt must also match configs/b2_training.json."
    )


if __name__ == "__main__":
    main()
