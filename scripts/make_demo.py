#!/usr/bin/env python3
"""Deterministic synthetic smoke-test inputs; none of these arrays are biological evidence."""
import json
from pathlib import Path
import numpy as np
import pandas as pd


def design(output, seed=42):
    output = Path(output)
    rng = np.random.default_rng(seed)
    output.mkdir(parents=True, exist_ok=True)
    x = rng.normal(size=(48, 16))
    rotation, _ = np.linalg.qr(rng.normal(size=(16, 16)))
    y = x @ rotation + rng.normal(scale=0.1, size=x.shape)
    query = rng.normal(size=(10, 16))
    # A nearby pair exercises the family/residual branch, rather than singleton-only output.
    query[1] = query[0] + rng.normal(scale=0.005, size=16)
    pool = rng.normal(size=(160, 16))
    projection = rng.normal(size=(16, 8))
    np.savez_compressed(
        output / "synthetic_inputs.npz",
        train_human=x,
        train_yeast=y,
        queries=query,
        yeast_pool=pool,
        synthetic_projection=projection,
    )
    # This small ridge fit demonstrates the design step; it is not the historical B2 training.
    w = np.linalg.solve(x.T @ x + np.eye(16), x.T @ y)
    np.savez_compressed(output / "demo_alignment.npz", ridge_weight=w)

    def norm(a):
        return a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)

    def pc(a):
        centered = a - a.mean(0)
        _, _, v = np.linalg.svd(centered, full_matrices=False)
        return centered - (centered @ v[:1].T) @ v[:1]

    members = [
        norm(pc(query)) @ norm(pc(pool)).T,
        norm(query @ w) @ norm(pool).T,
        norm(query @ projection) @ norm(pool @ projection).T,
    ]
    scores = np.zeros((10, 160))
    for member in members:
        order = np.argsort(-member, axis=1, kind="stable")
        ranks = np.empty_like(order)
        np.put_along_axis(ranks, order, np.arange(1, 161)[None, :], axis=1)
        scores += 1 / (20 + ranks)
    tasks = output / "tasks_ej"
    tasks.mkdir()
    genes = np.array([f"DEMO_Y{i:03d}" for i in range(160)])
    for i, row in enumerate(scores):
        order = np.argsort(-row, kind="stable")
        pd.DataFrame(
            {
                "rank": np.arange(1, 161),
                "yeast_gene": genes[order],
                "ensemble_score": row[order],
            }
        ).to_csv(tasks / f"yeast_task_DEMO_TARGET_{i:02d}.tsv", sep="\t", index=False)
    (tasks / "export_record.json").write_text(
        json.dumps(
            {
                "mode": "synthetic-demo",
                "seed": seed,
                "pool": 160,
                "tail_strains": [],
                "n_targets": 10,
                "disclosure": "synthetic embeddings; no ESM-2 or B2 model inference",
            }
        )
    )
    # Synthetic positive controls exercise the family/residual branches; they are not biological observations.
    compounds = pd.DataFrame(
        {
            "inchikey": ["DEMO_COMPOUND_A", "DEMO_COMPOUND_B", "DEMO_COMPOUND_C"],
            "smiles": ["CCO", "CC(=O)O", "C1=CC=CC=C1"],
            "pubchem_cid": ["", "", ""],
        }
    )
    compounds.to_csv(output / "compounds.tsv", sep="\t", index=False)
    centered = scores - scores.mean(axis=0)
    u, values, vt = np.linalg.svd(centered, full_matrices=False)
    dc = centered - np.outer(u[:, 0] * values[0], vt[0])
    response = rng.normal(size=(160, 6)).astype(np.float32)
    response[:, 0] = (dc[0] / np.std(dc[0]) + rng.normal(scale=0.05, size=160)).astype(
        np.float32
    )
    np.savez_compressed(
        output / "response.npz",
        strain_orfs=genes,
        compound_inchikeys=np.repeat(compounds.inchikey.to_numpy(dtype=str), 2),
        doses=np.tile(np.array(["1", "2"]), 3),
        is_vehicle=np.zeros(6, dtype=bool),
        z_score=response,
        method=np.array(["synthetic-demo-only"]),
    )
    return tasks, output / "response.npz", output / "compounds.tsv"
