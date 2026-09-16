#!/usr/bin/env python3
"""V6 neighborhood derivation (frozen inputs only; no labels/splits touched).

Neighborhood = descendant closure (is_a / part_of) of the pre-declared seeds
GO:0023052 (signaling) and GO:0055085 (transmembrane transport) in the frozen
go-basic.obo (go_bridge_assets_v1/raw), intersected with the V5.1-frozen
shared GO-BP ES term set, recomputed here with the identical recipe
(gene_bp_ancestors.tsv.gz + K562 73-column gate + es_term_min/max_genes 10/500
+ yeast matrix gene scope). Asserts the recomputation yields exactly 639
shared terms before intersecting.

Outputs (v6_process_bridge_20260912/frozen/):
  neighborhood_terms.tsv   term_id, name, seeds, n_human_genes, n_yeast_genes
  derivation_record.json   seeds, obo release/sha, input shas, counts, script sha
"""
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STAGE = Path("/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901")
MVP = STAGE / "yeastbridge_re_mvp"
OUT = STAGE / "v6_process_bridge_20260912/frozen"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MVP / "scripts"))
CFG = json.loads((MVP / "configs/v5_1_joint_formal.json").read_text())
SEEDS = {"GO:0023052": "signaling", "GO:0055085": "transmembrane_transport"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ---- 1) exact V5.1 shared-term recomputation ----
import anndata as ad  # noqa: E402

gate = CFG["data_quality_gate"]["column_abs_max_threshold"]
a = ad.read_h5ad(STAGE / CFG["inputs"]["k562_h5ad"]["path"])
genes_k = list(a.var.gene_name)
obs = list(a.obs.index)
pat = re.compile(r"^\d+_([^_]+)_(P1|P1P2)_")
keep = [n for n in obs if pat.match(n) and pat.match(n).group(2) in {"P1", "P2", "P1P2"}]
Xp = a.X[[obs.index(n) for n in keep]]
Xp = Xp.toarray() if hasattr(Xp, "toarray") else np.asarray(Xp)
Xp = np.nan_to_num(Xp, nan=0.0)
col_ok = np.abs(Xp).max(0) < gate
n_drop = int((~col_ok).sum())
assert n_drop == CFG["data_quality_gate"]["expected_dropped_columns"], n_drop
genes_k = [g for g, ok in zip(genes_k, col_ok) if ok]
gidx_k = set(genes_k)
del a, Xp

from cross_species_match import load_matrix  # noqa: E402

genes_y, mat_y = load_matrix(STAGE / CFG["inputs"]["yeast_matrix"]["path"])
gidx_y = set(genes_y)
del mat_y

anc = pd.read_csv(STAGE / CFG["inputs"]["go_bp_ancestors"]["path"], sep="\t",
                  compression="gzip", usecols=["species", "gene", "ancestor_go_id"])
tmin = CFG["representation"]["es_term_min_genes"]
tmax = CFG["representation"]["es_term_max_genes"]


def term_sets(df):
    d = {}
    for g, t in zip(df.gene, df.ancestor_go_id):
        d.setdefault(t, []).append(g)
    return {t: gs for t, gs in d.items() if tmin <= len(gs) <= tmax}


ts_h = term_sets(anc[(anc.species == "human") & anc.gene.isin(gidx_k)])
ts_y = term_sets(anc[(anc.species == "yeast") & anc.gene.isin(gidx_y)])
shared = sorted(set(ts_h) & set(ts_y))
print("shared terms:", len(shared), flush=True)
assert len(shared) == 639, f"expected 639 shared terms, got {len(shared)}"

# ---- 2) obo descendant closure of seeds ----
obo = STAGE / "go_bridge_assets_v1/raw/go-basic.obo"
name_of = {}
parents = {}
cur_id = None
in_term = False
release = None
for line in obo.open():
    line = line.rstrip("\n")
    if line.startswith("data-version:") and release is None:
        release = line.split(":", 1)[1].strip()
    if line == "[Term]":
        in_term = True
        cur_id = None
        continue
    if line.startswith("[") and line != "[Term]":
        in_term = False
        continue
    if not in_term:
        continue
    if line.startswith("id: "):
        cur_id = line[4:].strip()
        parents.setdefault(cur_id, set())
    elif cur_id is None:
        continue
    elif line.startswith("name: "):
        name_of[cur_id] = line[6:].strip()
    elif line.startswith("is_a: "):
        parents[cur_id].add(line[6:].split()[0].strip())
    elif line.startswith("relationship: part_of "):
        parents[cur_id].add(line[len("relationship: part_of "):].split()[0].strip())

children = {}
for t, ps in parents.items():
    for p in ps:
        children.setdefault(p, set()).add(t)

desc = {}
for s in SEEDS:
    assert s in parents, f"seed {s} not found in obo"
    stack, seen = [s], {s}
    while stack:
        t = stack.pop()
        for c in children.get(t, ()):
            if c not in seen:
                seen.add(c)
                stack.append(c)
    desc[s] = seen
    print(f"seed {s} ({name_of.get(s)}): {len(seen)} descendants", flush=True)

nbr = sorted((desc[SEEDS and "GO:0023052"] | desc["GO:0055085"]) & set(shared))
print("neighborhood terms (∩ shared):", len(nbr), flush=True)

rows = []
for t in nbr:
    seeds = "+".join(nm for s, nm in SEEDS.items() if t in desc[s])
    rows.append({"term_id": t, "name": name_of.get(t, ""), "seeds": seeds,
                 "n_human_genes": len(ts_h[t]), "n_yeast_genes": len(ts_y[t])})
pd.DataFrame(rows).to_csv(OUT / "neighborhood_terms.tsv", sep="\t", index=False)

rec = {"seeds": SEEDS, "obo_path": str(obo), "obo_sha256": sha(obo),
       "obo_release": release, "shared_terms": len(shared),
       "neighborhood_terms": len(nbr),
       "inputs_sha256": {k: sha(STAGE / v["path"]) for k, v in CFG["inputs"].items()
                         if k in ("k562_h5ad", "yeast_matrix", "go_bp_ancestors")},
       "script_sha256": sha(Path(__file__).resolve())}
(OUT / "derivation_record.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))
print(json.dumps({k: rec[k] for k in ("obo_release", "shared_terms", "neighborhood_terms")},
                 ensure_ascii=False))
print("[done]", OUT)
