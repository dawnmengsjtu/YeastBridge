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


def check_destination(target, expected):
    if target.exists() and sha256(target) != expected:
        raise ValueError(f"Refusing to replace a different local asset: {target}")


def install_verified(source, target, expected):
    check_destination(target, expected)
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    pending = target.with_name(target.name + ".installing")
    try:
        shutil.copyfile(source, pending)
        if sha256(pending) != expected:
            raise ValueError(f"Copied asset checksum mismatch: {target}")
        pending.replace(target)
    finally:
        pending.unlink(missing_ok=True)


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
    cells_source = ROOT / "models/training/evidence/corpus_cells.tsv"
    cells_sha = cfg["required_assets"]["cells"]["sha256"]
    if sha256(cells_source) != cells_sha:
        p.error("Cell roster checksum mismatch")
    # Check every existing destination before installing any asset.
    for name, expected_sha in expected_files.items():
        check_destination(destination / name, expected_sha)
    check_destination(destination / "corpus_cells.tsv", cells_sha)
    if args.checkpoint:
        expected_checkpoint = cfg["required_assets"]["pretrained_checkpoint"]["sha256"]
        if sha256(args.checkpoint) != expected_checkpoint:
            p.error("Pretrained checkpoint checksum mismatch; no model was installed")
        checkpoint_target = ROOT / cfg["required_assets"]["pretrained_checkpoint"]["path"]
        check_destination(checkpoint_target, expected_checkpoint)
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
            install_verified(Path(temp) / name, destination / name, expected_files[name])
    install_verified(cells_source, destination / "corpus_cells.tsv", cells_sha)
    if args.checkpoint:
        install_verified(args.checkpoint, checkpoint_target, expected_checkpoint)
    print(
        "Training inputs verified. The upstream models.ckpt must also match configs/b2_training.json."
    )


if __name__ == "__main__":
    main()
