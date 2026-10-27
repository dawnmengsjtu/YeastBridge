#!/usr/bin/env python3
"""Replay the registered SGA graph-native T2 and T4 comparator protocols."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import paper_release as paper
from scripts.model_selection.step4_route_confirmation.eval import data as ed
ALPHA=0.3

def _load_sga_neighbors(root: Path, path: str) -> dict[str, set[str]]:
    file = root / path
    neighbors: dict[str, set[str]] = {}
    import gzip

    with gzip.open(file, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query, array = str(row.get("query_orf", "")), str(row.get("array_orf", ""))
            if query and array:
                neighbors.setdefault(query, set()).add(array)
                neighbors.setdefault(array, set()).add(query)
    return neighbors

def _propagate_personalized_pagerank(
    nodes: list[str],
    neighbors: dict[str, set[str]],
    anchor_weights: dict[str, float],
    alpha: float = 0.30,
    tolerance: float = 1e-10,
    max_iterations: int = 300,
) -> tuple[np.ndarray, int]:
    """Signed personalised PageRank over the undirected SGA graph.

    Anchor mass restarts every iteration (p = alpha*anchor + (1-alpha)*W^T p),
    so unlike repeated neighbour-mean diffusion the signal cannot wash out to
    uniformity on a dense graph. Positive and negative anchor masses propagate
    as separate channels and recombine with their signs, preserving
    intervention direction.
    """

    import scipy.sparse as sp

    index_of = {node: i for i, node in enumerate(nodes)}
    size = len(nodes)
    rows, cols = [], []
    seen_pairs: set[tuple[int, int]] = set()
    for orf, connected in neighbors.items():
        i = index_of.get(orf)
        if i is None:
            continue
        for other in connected:
            j = index_of.get(other)
            if j is None or (i, j) in seen_pairs:
                continue
            seen_pairs.add((i, j))
            rows.append(i)
            cols.append(j)
    data = np.ones(len(rows), dtype=np.float64)
    adjacency = sp.coo_matrix((data, (rows, cols)), shape=(size, size)).tocsr()
    row_sums = np.asarray(adjacency.sum(axis=1)).ravel()
    walk = sp.diags(1.0 / np.maximum(row_sums, 1.0)) @ adjacency

    def run_channel(mass: np.ndarray) -> tuple[np.ndarray, int]:
        restart = alpha * mass
        value = restart.copy()
        iterations = max_iterations
        for step in range(1, max_iterations + 1):
            nxt = restart + (1.0 - alpha) * (walk.T @ value)
            delta = float(np.abs(nxt - value).max())
            value = nxt
            if delta < tolerance:
                iterations = step
                break
        return value, iterations

    def channel_vector(channel: dict[str, float]) -> np.ndarray:
        mass = np.zeros(size)
        total = sum(abs(w) for w in channel.values())
        if total <= 0.0:
            return mass
        for orf, w in channel.items():
            i = index_of.get(orf)
            if i is not None:
                mass[i] = abs(w) / total
        return mass

    positives = {o: w for o, w in anchor_weights.items() if w > 0}
    negatives = {o: w for o, w in anchor_weights.items() if w < 0}
    pos_score, iter_pos = run_channel(channel_vector(positives))
    neg_score, iter_neg = run_channel(channel_vector(negatives))

    pos_scale = (
        sum(w for w in positives.values()) / len(positives) if positives else 0.0
    )
    neg_scale = (
        abs(sum(w for w in negatives.values())) / len(negatives) if negatives else 0.0
    )
    return pos_scale * pos_score - neg_scale * neg_score, max(iter_pos, iter_neg)

def propagate(seeds, nodes, neighbors, alpha):
    vec, _ = _propagate_personalized_pagerank(nodes, neighbors, seeds, alpha=alpha)
    return vec

def e_t2_essentiality(nodes, index_of, neighbors) -> dict:
    lab = ed.essentiality()
    lab = lab[lab["essentiality"].isin(["essential", "nonessential"])]
    ess = [g for g in lab.loc[lab["essentiality"] == "essential", "systematic"] if g in index_of]
    non = [g for g in lab.loc[lab["essentiality"] == "nonessential", "systematic"] if g in index_of]
    print(f"[E-T2] essential on graph={len(ess)} nonessential pool={len(non)}", flush=True)
    non_idx = np.array([index_of[g] for g in non])
    ranks = []
    for k, o in enumerate(ess):
        seeds = {g: 1.0 for g in ess if g != o}
        vec = propagate(seeds, nodes, neighbors, ALPHA)
        s_pool = vec[non_idx]
        s_o = vec[index_of[o]]
        r = 1.0 + np.count_nonzero(s_pool > s_o) + np.count_nonzero(s_pool == s_o) / 2.0
        ranks.append(r / len(non))
        if (k + 1) % 200 == 0:
            print(f"[E-T2] {k+1}/{len(ess)} median={np.median(ranks):.4f}", flush=True)
    ranks = np.array(ranks)
    return {"median_normalized_rank": round(float(np.median(ranks)), 5),
            "pooled_auroc": round(float(np.mean(1.0 - ranks)), 5),
            "n_essential_loo": len(ess), "n_pool": len(non)}

def e_t4_engineering(nodes, index_of, neighbors) -> dict:
    rec = ed.engineering_records()
    pw = ed.pathway_genes()
    hits3 = hits5 = dir_ok = n = 0
    for _, r in rec.iterrows():
        targets = [g.strip() for g in str(r["systematic"]).split("/") if g.strip()]
        pid = r["pathway_id"]
        pw_genes = pw.loc[pw["pathway_id"] == pid, "systematic"].tolist()
        universe = sorted(set(pw_genes) | set(targets))
        if r["direction"] not in ("overexpress", "knockdown"):
            continue
        for g in targets:
            seeds = {m: 1.0 for m in universe if m != g and m in index_of}
            if not seeds or g not in index_of:
                continue
            vec = propagate(seeds, nodes, neighbors, ALPHA)
            sc = {m: float(vec[index_of[m]]) for m in universe if m in index_of}
            order_oe = sorted(sc, key=lambda x: sc[x])
            order_kd = sorted(sc, key=lambda x: -sc[x])
            order = order_oe if r["direction"] == "overexpress" else order_kd
            rank = order.index(g) + 1
            n += 1
            hits3 += rank <= 3
            hits5 += rank <= 5
            dir_ok += (sc[g] < 0) == (r["direction"] == "overexpress")
    return {"n_usable": n, "hit_at_3": round(hits3 / max(n, 1), 4),
            "hit_at_5": round(hits5 / max(n, 1), 4),
            "direction_consistency": round(dir_ok / max(n, 1), 4)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    paper.require('graph');paper.require('comparisons')
    output=paper.new_output(a.output)
    neighbors=_load_sga_neighbors(ROOT,'raw/model_selection/graph/sga.corrected.tsv.gz')
    nodes=sorted({n for g,connected in neighbors.items() for n in [g,*connected]})
    index={n:i for i,n in enumerate(nodes)}
    result=dict(T2=e_t2_essentiality(nodes,index,neighbors),T4=e_t4_engineering(nodes,index,neighbors))
    paper.write_json(output/'result.json',result)
    reference=json.loads((ROOT/'raw/model_selection/reference/step4_route_confirmation/five_route_results.json').read_text())['tasks']['E_sga_propagation']
    checks={task:result[task]==reference[task] for task in ['T2','T4']}
    paper.write_json(output/'verification.json',dict(status='passed' if all(checks.values()) else 'different',graph_nodes=len(nodes),metrics=result,checks=checks,protocol='Original graph-native T2/T4; no embedding substitution',script_sha256=paper.digest(Path(__file__))))
    print(json.dumps(result,indent=2))
    return 0 if all(checks.values()) else 2

if __name__=='__main__':raise SystemExit(main())
