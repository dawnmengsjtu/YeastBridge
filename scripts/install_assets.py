#!/usr/bin/env python3
"""Download or install large assets after frozen checksum verification."""
import argparse
import json
import shutil
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256

ASSET_MANIFEST = json.loads((ROOT / "data/asset_downloads.json").read_text())["assets"]
ASSETS = {
    key: (ASSET_MANIFEST[key]["path"], ASSET_MANIFEST[key]["sha256"])
    for key in ("weights", "response")
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights", type=Path)
    ap.add_argument("--response", type=Path)
    ap.add_argument(
        "--download", action="store_true", help="Download the two frozen runtime assets"
    )
    ap.add_argument(
        "--legacy",
        action="store_true",
        help="Install the archived September B2 checkpoint",
    )
    args = ap.parse_args()
    weight_key = (
        "weights_legacy"
        if args.legacy and "weights_legacy" in ASSET_MANIFEST
        else "weights"
    )
    assets = dict(ASSETS)
    assets["weights"] = (
        ASSET_MANIFEST[weight_key]["path"],
        ASSET_MANIFEST[weight_key]["sha256"],
    )
    if args.download:
        if args.weights or args.response:
            ap.error("--download cannot be combined with local file arguments")
        from scripts.download_asset import download

        for key in assets:
            download(weight_key if key == "weights" else key)
        return
    supplied = [(k, getattr(args, k)) for k in assets if getattr(args, k) is not None]
    if not supplied:
        ap.error("supply --download, --weights and/or --response")
    # Validate every source before copying any of them.
    for key, source in supplied:
        if sha256(source) != assets[key][1]:
            ap.error(f"{key}: checksum mismatch; refusing to install a different asset")
    for key, source in supplied:
        destination = ROOT / assets[key][0]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.resolve() != destination.resolve():
            shutil.copyfile(source, destination)
        print(f"Installed and verified: {assets[key][0]}")


if __name__ == "__main__":
    main()
