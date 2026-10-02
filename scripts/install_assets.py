#!/usr/bin/env python3
"""Download or install large assets after frozen checksum verification."""
import argparse
import shutil
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256

ASSETS = {
    "weights": (
        "raw/tier1_models/route_b/final_model.pt",
        "112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7",
    ),
    "response": (
        "raw/tier1_response/strain_response.npz",
        "d174dc0a1cc403b9ab4ce730aa07a11eb51aa290fe04196c2c903d47c486aadf",
    ),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights", type=Path)
    ap.add_argument("--response", type=Path)
    ap.add_argument(
        "--download", action="store_true", help="Download the two frozen runtime assets"
    )
    args = ap.parse_args()
    if args.download:
        if args.weights or args.response:
            ap.error("--download cannot be combined with local file arguments")
        from scripts.download_asset import download

        for key in ASSETS:
            download(key)
        return
    supplied = [(k, getattr(args, k)) for k in ASSETS if getattr(args, k) is not None]
    if not supplied:
        ap.error("supply --download, --weights and/or --response")
    # Validate every source before copying any of them.
    for key, source in supplied:
        if sha256(source) != ASSETS[key][1]:
            ap.error(f"{key}: checksum mismatch; refusing to install a different asset")
    for key, source in supplied:
        destination = ROOT / ASSETS[key][0]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.resolve() != destination.resolve():
            shutil.copyfile(source, destination)
        print(f"Installed and verified: {ASSETS[key][0]}")


if __name__ == "__main__":
    main()
