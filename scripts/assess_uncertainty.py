#!/usr/bin/env python3
"""Report conditional Monte Carlo uncertainty of candidate permutation p values."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    sources = {}
    for direction in ["pos", "neg"]:
        matrix = args.run_dir / f"confirm_{direction}/exec_matrix.tsv"
        summary = matrix.with_name("execute_summary.json")
        data = pd.read_csv(matrix, sep="\t")
        sources[str(matrix)] = sha256(matrix)
        sources[str(summary)] = sha256(summary)
        if data.empty:
            continue
        n = int(json.loads(summary.read_text())["n_perm"])
        for row in data.itertuples():
            b_float = float(row.emp_p) * (n + 1) - 1
            b = round(b_float)
            if not 0 <= b <= n or not np.isclose(b_float, b, atol=1e-7, rtol=0):
                raise ValueError(
                    "p value does not match the recorded permutation count"
                )
            interval = binomtest(b, n).proportion_ci(
                confidence_level=0.95, method="exact"
            )
            rows.append(
                {
                    "target_id": row.target_id,
                    "inchikey": row.inchikey,
                    "direction": direction,
                    "emp_p": float(row.emp_p),
                    "n_perm": n,
                    "null_exceedances": b,
                    "null_probability_ci95_low": interval.low,
                    "null_probability_ci95_high": interval.high,
                }
            )
    out = args.run_dir / "uncertainty"
    out.mkdir(exist_ok=False)
    table = out / "monte_carlo_intervals.tsv"
    pd.DataFrame(rows).to_csv(table, sep="\t", index=False)
    write_json(
        out / "interpretation.json",
        {
            "rows": len(rows),
            "confidence_level": 0.95,
            "method": "Clopper-Pearson exact binomial interval",
            "quantity": "Probability that a random strain permutation produces a maximum-over-dose statistic at least as large as observed",
            "emp_p_formula": "(null_exceedances + 1)/(n_perm + 1)",
            "interval_formula": "Exact binomial confidence interval for b of B exceedances, before the +1 p-value correction",
            "scope": "Monte Carlo sampling only, conditional on the selected candidate and observed assay data",
            "not_covered": [
                "stage-one nomination selection",
                "assay measurement uncertainty",
                "biological target specificity",
                "drug potency",
            ],
            "source_sha256": sources,
            "output_sha256": sha256(table),
        },
    )
    print(f"Conditional Monte Carlo intervals: {len(rows)} candidates")


if __name__ == "__main__":
    main()
