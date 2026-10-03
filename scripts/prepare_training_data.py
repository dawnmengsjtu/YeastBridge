#!/usr/bin/env python3
"""Rebuild B2 inputs from the original GEO count TSV and frozen ESM-2 embeddings."""
import argparse
import json
from pathlib import Path
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts-tsv",
        type=Path,
        required=True,
        help="Decompressed GEO GSE125162_ALL-fastqTomat0-Counts.tsv",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="New or empty output directory"
    )
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("Output must be a new or empty directory")
    recovery = json.loads((ROOT / "models/training/evidence/recovery.json").read_text())
    if sha256(args.counts_tsv) != recovery["raw_counts"]["sha256"]:
        parser.error("Raw GEO count TSV differs from the recorded snapshot")
    master_path = ROOT / "raw/externalvalidation/mappings/gene_master.tsv"
    # Git's historical LF copy and the original CRLF file contain identical rows.
    master_sha = sha256(master_path)
    if master_sha not in {
        recovery["gene_master"]["sha256"],
        "985035e0ed520eca266f69e343d05e0d0b84dc87a7fa4a620623beb5c73dfa07",
    }:
        parser.error("Gene master differs from the recorded snapshot")
    master = pd.read_csv(master_path, sep="\t", dtype=str).fillna("")
    genes = [g.strip() for g in master.systematic if g.strip()]
    esm_path = ROOT / "raw/tier1_esm2/yeast_650m/esm2_mean_fp32.npy"
    index_path = esm_path.with_name("index.tsv")
    frozen = json.loads((ROOT / "configs/v7_b2_joint_formal.json").read_text())[
        "inputs"
    ]["yeast_esm2"]
    if (
        sha256(esm_path) != frozen["sha256"]
        or sha256(index_path) != frozen["index_sha256"]
    ):
        parser.error(
            "Yeast protein embedding or index differs from the frozen snapshot"
        )
    index = pd.read_csv(index_path, sep="\t", dtype=str).fillna("")
    rows = {g: i for i, g in enumerate(index.systematic)}
    esm = np.load(esm_path, mmap_mode="r")
    protein = np.zeros((len(genes), esm.shape[1]), dtype=np.float32)
    for i, g in enumerate(genes):
        if g in rows:
            protein[i] = esm[rows[g]]
    args.output.mkdir(parents=True, exist_ok=True)
    np.save(args.output / "B2_esm2_matrix.npy", protein)
    # A streaming pass bounds memory without changing the original row order.
    with args.counts_tsv.open() as file:
        n_cells = sum(1 for _ in file) - 1
    matrix = np.lib.format.open_memmap(
        args.output / "corpus_counts.npy",
        mode="w+",
        dtype=np.float32,
        shape=(n_cells, len(genes)),
    )
    matrix[:] = 0
    cells = []
    offset = 0
    hit = None
    for frame in pd.read_csv(args.counts_tsv, sep="\t", index_col=0, chunksize=2048):
        columns = {c: i for i, c in enumerate(frame.columns)}
        matched = [(j, columns[g]) for j, g in enumerate(genes) if g in columns]
        hit = len(matched)
        for j, k in matched:
            matrix[offset : offset + len(frame), j] = frame.iloc[:, k].to_numpy(
                dtype=np.float32
            )
        cells.extend(frame.index.tolist())
        offset += len(frame)
    if offset != n_cells:
        raise ValueError("TSV row count mismatch")
    matrix.flush()
    del matrix
    pd.DataFrame({"cell": cells}).to_csv(
        args.output / "corpus_cells.tsv", sep="\t", index=False
    )
    cfg = json.loads((ROOT / "configs/b2_training_v2.json").read_text())
    checks = {}
    for name in ["counts", "protein_matrix", "cells"]:
        item = cfg["required_assets"][name]
        path = args.output / Path(item["path"]).name
        checks[path.name] = {"sha256": sha256(path), "expected_sha256": item["sha256"]}
        checks[path.name]["match"] = checks[path.name]["sha256"] == item["sha256"]
    record = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": sha256(args.counts_tsv),
        "gene_master_sha256": master_sha,
        "script_sha256": sha256(Path(__file__)),
        "n_cells": n_cells,
        "n_genes": len(genes),
        "count_columns_matched": hit,
        "missing_protein_rows": sum(g not in rows for g in genes),
        "checks": checks,
        "status": "pass" if all(c["match"] for c in checks.values()) else "fail",
    }
    write_json(args.output / "preprocessing.json", record)
    print(json.dumps(record, indent=2))
    return 0 if record["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
