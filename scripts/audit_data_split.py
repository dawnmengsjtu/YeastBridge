#!/usr/bin/env python3
"""Audit the frozen supervision partition and optional B2 count-row duplicates."""
import argparse
import hashlib
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path)
    parser.add_argument(
        "--cells", type=Path, default=ROOT / "models/training/evidence/corpus_cells.tsv"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = ROOT / "results_a6000/v51_per_query_ranks_trackA.tsv"
    pairs = pd.read_csv(source, sep="\t")
    overlaps = {}
    for key in ["human", "yeast_orf"]:
        groups = {label: set(frame[key]) for label, frame in pairs.groupby("split")}
        overlaps[key] = {
            f"{a}|{b}": sorted(groups[a] & groups[b])
            for a, b in combinations(groups, 2)
        }
    record = {
        "supervision": {
            "sha256": sha256(source),
            "rows": len(pairs),
            "splits": pairs.groupby("split").size().to_dict(),
            "duplicate_pairs": int(pairs.duplicated(["human", "yeast_orf"]).sum()),
            "cross_partition_overlap": overlaps,
        },
        "scope": "Exact identifiers and exact count-row duplication; not a claim of biological or hidden-test independence",
    }
    if args.counts:
        counts = np.load(args.counts, mmap_mode="r")
        cells = pd.read_csv(args.cells, sep="\t")
        if len(cells) != len(counts):
            raise ValueError("Cell IDs and count rows disagree")
        order = np.random.default_rng(42).permutation(len(counts))
        validation = set(order[: int(0.05 * len(counts))].tolist())
        hashes = {}
        duplicates = 0
        cross = 0
        for i in range(len(counts)):
            digest = hashlib.sha256(counts[i].tobytes()).hexdigest()
            split = "validation" if i in validation else "train"
            if digest in hashes:
                duplicates += 1
                cross += int(hashes[digest] != split)
            else:
                hashes[digest] = split
        record["b2"] = {
            "counts_sha256": sha256(args.counts),
            "cells_sha256": sha256(args.cells),
            "shape": list(counts.shape),
            "duplicate_cell_ids": int(cells.cell.duplicated().sum()),
            "duplicate_count_rows": duplicates,
            "cross_partition_duplicate_count_rows": cross,
            "seed": 42,
            "train_cells": len(counts) - len(validation),
            "validation_cells": len(validation),
            "split_unit": "cell",
            "batch_grouped": False,
            "independent_test_set": False,
        }
    problems = [
        items for values in overlaps.values() for items in values.values() if items
    ]
    record["status"] = (
        "pass"
        if not problems
        and not record["supervision"]["duplicate_pairs"]
        and not record.get("b2", {}).get("cross_partition_duplicate_count_rows", 0)
        else "fail"
    )
    write_json(args.output, record)
    print(f"{record['status']}: {args.output}")
    return 0 if record["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
