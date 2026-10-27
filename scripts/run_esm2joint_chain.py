#!/usr/bin/env python3
"""Compatibility entry for the manuscript's frozen-model reconstruction.

Usage: python scripts/run_esm2joint_chain.py --output outputs/full
Requires both inference and numerical archives. Rebuilds profiles, verifies
them against the screening inputs, and reconstructs the screening decisions.
The original competition dispatcher is preserved in its historical release.
"""
import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("outputs/full"))
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    if args.threads < 1:
        p.error("--threads must be positive")
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = str(args.threads)
    from scripts import paper_release as paper
    paper.require("all")
    output = paper.new_output(args.output.resolve())
    report = {
        "profiles": paper.profiles(output / "profiles", args.threads),
        "numerical": paper.reproduce(output / "numerical"),
    }
    paper.write_json(output / "verification.json", report)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2)
