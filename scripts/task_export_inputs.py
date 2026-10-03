"""Load only the representations and training pairs used by EJ task export.

The historical v7_setup also loads scYeast, expression benchmarks and GO labels
for its separate evaluation arms. None of those enters an EJ candidate score.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from project_io import ROOT, sha256
from v5_joint_formal_run import l2n


def setup(config_path):
    cfg = json.loads(Path(config_path).read_text())
    root = (ROOT / cfg.get("stage_root", ".")).resolve()
    inputs = cfg["inputs"]
    for label, item in inputs.items():
        for pk, hk in [
            ("path", "sha256"),
            ("index_path", "index_sha256"),
            ("path_abs", "sha256"),
            ("index_path_abs", "index_sha256"),
        ]:
            if pk in item:
                if not item.get(hk) or sha256(root / item[pk]) != item[hk]:
                    raise ValueError(f"Input checksum mismatch: {label}.{pk}")

    def path(label, key="path"):
        return root / inputs[label][key]

    human = np.load(path("esm2_human")).astype(np.float64)
    human_index = pd.read_csv(path("esm2_human", "index_path"), sep="\t").fillna("")
    sym2row = {}
    for i, symbol in enumerate(human_index.common.astype(str)):
        sym2row.setdefault(symbol, i)
    checkpoint = torch.load(
        path("route_b_model"), map_location="cpu", weights_only=False
    )
    state = checkpoint["state_dict"]
    weight = state["pos_emb.proj.weight"].float().numpy().astype(np.float64)
    bias = state["pos_emb.proj.bias"].float().numpy().astype(np.float64)
    table = np.load(path("route_b_table")).astype(np.float64)
    order = pd.read_csv(
        path("route_b_gene_order"), sep="\t", dtype=str
    ).systematic.tolist()
    rb_row = {gene: i for i, gene in enumerate(order)}
    yeast = np.load(path("yeast_esm2", "path_abs")).astype(np.float64)
    yeast_index = pd.read_csv(path("yeast_esm2", "index_path_abs"), sep="\t")
    yeast_row = {gene: i for i, gene in enumerate(yeast_index.systematic.astype(str))}
    metadata = pd.read_csv(path("yeast_metadata"), sep="\t")
    metadata = metadata[
        metadata.benchmark_eligible.astype(str).str.lower() == "1"
    ].reset_index(drop=True)
    pool = [
        str(g)
        for g in metadata.target_stable_id
        if str(g) in rb_row and str(g) in yeast_row
    ]
    pool_pos = {gene: i for i, gene in enumerate(pool)}
    pairs = pd.read_csv(path("v51_pair_list"), sep="\t")
    if len(pairs) != 753:
        raise ValueError("Expected the frozen 753-pair development partition")
    prim = [
        (str(r.human), str(r.yeast_orf), str(r.split), int(r.fold))
        for r in pairs.itertuples()
        if str(r.human) in sym2row and str(r.yeast_orf) in pool_pos
    ]
    query = sorted({h for h, _, _, _ in prim})
    return {
        "X_es_all": l2n(human),
        "Y_es": l2n(yeast[[yeast_row[g] for g in pool]]),
        "Y_rb": l2n(table[[rb_row[g] for g in pool]]),
        "W_inj": weight,
        "b_inj": bias,
        "sym2row_h": sym2row,
        "meta_y": metadata,
        "pool": pool,
        "pool_pos": pool_pos,
        "prim": prim,
        "query_genes": query,
        "q_pos": {gene: i for i, gene in enumerate(query)},
    }
