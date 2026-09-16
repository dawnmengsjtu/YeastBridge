#!/usr/bin/env python3
"""V7-B2 formal asset build (2026-09-12, user-directed).

Human protein sequences for every gene the V7-B2 three-arm two-gate formal run
needs: 722 OrthoDB-pair genes (seed cache copied verbatim from
esm2_joint_test_20260912/build/uniprot_cache, 721 ok / 1 no_result NEDD8-MDP1)
+ 94 Norman perturbation genes (89 fetched here, resume-safe).

Fetcher conventions identical to esm2_joint_test build: UniProt REST
(reviewed, human, gene_exact), 3 retries with backoff, 0.12s polite delay.
Outputs: uniprot_cache/<SYM>.json, manifest.tsv, batch.fasta
(fasta header mirrors esm2_cross_match.py: >sp|ACC|SYM_HUMAN GN=SYM).
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

STAGE = Path("/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901")
ASSETS = STAGE / "v7_b2_formal_20260912/assets"
CACHE = ASSETS / "uniprot_cache"
CACHE.mkdir(parents=True, exist_ok=True)

edges = STAGE / "k562_bridge_input_v1/official_run/orthodb_bipartite_edges.tsv"
pair_genes = sorted({l.split("\t")[0] for l in edges.read_text().splitlines()[1:] if l.strip()})
mn = pd.read_csv(STAGE / "gse_go_three_arm_v3_open_set/inputs_v3_open_set/higher_metadata.tsv", sep="\t")
mn = mn[mn.benchmark_eligible.astype(str).str.lower() == "true"]
norman_genes = set()
for p in mn.perturbation_genes:
    norman_genes.update(x.strip() for x in re.split(r"[+;,|]", str(p)) if x.strip())
genes = sorted(set(pair_genes) | norman_genes)
print(f"genes total: {len(genes)} (pair {len(pair_genes)}, norman {len(norman_genes)})", flush=True)


def fetch(sym):
    q = urllib.parse.quote(f"gene_exact:{sym} AND organism_id:9606 AND reviewed:true")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession,gene_primary,sequence&format=json&size=1")
    for attempt, wait in enumerate((2, 5, 12)):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.loads(r.read())
            res = data.get("results", [])
            if not res:
                return {"symbol": sym, "status": "no_result"}
            hit = res[0]
            seq = hit.get("sequence", {}).get("value")
            acc = hit.get("primaryAccession")
            if not seq or not acc:
                return {"symbol": sym, "status": "no_sequence"}
            return {"symbol": sym, "status": "ok", "accession": acc,
                    "gene": (hit.get("genes") or [{}])[0].get("geneName", {}).get("value", sym),
                    "seq_len": len(seq), "sequence": seq}
        except Exception as exc:  # noqa: BLE001
            if attempt == 2:
                return {"symbol": sym, "status": f"error:{type(exc).__name__}"}
            time.sleep(wait)
    return {"symbol": sym, "status": "error:unreachable"}


n_fetch = 0
for sym in genes:
    fp = CACHE / f"{sym}.json"
    if fp.exists():
        continue
    d = fetch(sym)
    fp.write_text(json.dumps(d))
    n_fetch += 1
    time.sleep(0.12)
    if n_fetch % 10 == 0:
        print(f"fetched {n_fetch}...", flush=True)
print(f"newly fetched: {n_fetch}", flush=True)

pg, ng_ = set(pair_genes), norman_genes
rows = []
for sym in genes:
    fp = CACHE / f"{sym}.json"
    d = json.loads(fp.read_text()) if fp.exists() else {"symbol": sym, "status": "missing"}
    rows.append({"symbol": sym, "status": d.get("status"), "accession": d.get("accession", ""),
                 "seq_len": d.get("seq_len", ""), "in_pairs": int(sym in pg), "in_norman": int(sym in ng_)})
pd.DataFrame(rows).to_csv(ASSETS / "manifest.tsv", sep="\t", index=False)
n_ok = sum(r["status"] == "ok" for r in rows)
with (ASSETS / "batch.fasta").open("w") as fh:
    for r in rows:
        if r["status"] != "ok":
            continue
        d = json.loads((CACHE / (r["symbol"] + ".json")).read_text())
        fh.write(">sp|{}|{}_HUMAN GN={}\n{}\n".format(
            d["accession"], r["symbol"], r["symbol"], d["sequence"]))
print(f"manifest: ok={n_ok} bad={len(rows) - n_ok}; batch.fasta written", flush=True)
