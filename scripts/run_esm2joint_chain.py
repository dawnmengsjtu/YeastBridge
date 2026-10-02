#!/usr/bin/env python3
"""Compatibility entry for the portable complete EJ-dc pipeline.

Usage: python scripts/run_esm2joint_chain.py --output outputs/full
The original server-specific dispatcher is preserved in Git history.
"""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from main import pipeline

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("outputs/full"))
    args = p.parse_args()
    try:
        pipeline("full", args.output.resolve())
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2)
