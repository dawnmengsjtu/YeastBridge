# Real-data preflight — 2026-09-01

The server already contains two genuine signed perturbation sources:

- Norman CRISPRa: 223 human single/double perturbations × 3,066 response genes.
- SGA: 977,412 signed deletion/DAmP genetic-interaction edges covering 5,113
  query ORFs and 4,451 array ORFs.

The 977,412-row curated SGA file is a thresholded (`p < 0.05`,
`abs(score) >= 0.08`) sparse view. Its absent pairs are unknown rather than zero.
For modelling, use the four complete DAmP, ExE, EXN/NXE and NxN pairwise tables
(20,705,612 rows total); reserve the curated file for fast QC. All four expose
`score`, `pvalue`, query/array single-mutant fitness, double-mutant fitness and
its standard deviation.

They are scientifically eligible input types, but the current bridge is not ready
for a real three-arm claim. With the frozen OrthoDB pairs and 47 yeast network
modules, only 362/3,066 Norman response genes map (11.8%). Of 100 single
perturbations, only six map to one unambiguous module, distributed 1/3/1/1 across
four modules. That is insufficient for leakage-resistant train/test evaluation.

The existing Norman development/final-test split also shares 42 perturbation genes.
Four connected gene components cross the partition; the largest component contains
57 genes. A new component-grouped split is required before testing combinations.

Lee 2014 HIP/HOP is present as a 5,668 × 3,850 signed strain-by-drug response
matrix. It is valuable downstream, but it is not KO truth and is explicitly
excluded from the cross-species benchmark labels. `becker_epi.h5ad` is human CRC
observational single-cell data, not yeast deletion data.

Run the local manifest check:

```bash
python scripts/cross_species_real_preflight.py \
  --manifest configs/cross_species_match.real_assets.json
```

After placing the script on the server or mounting the assets, add
`--verify-files` to verify sizes and supplied SHA-256 hashes. The source files are
never modified.

Fastest scientifically defensible unblock:

1. Freeze human and yeast GO Biological Process annotations plus a GO-slim or
   ontology-ancestor projection; keep OrthoDB edges at highest confidence.
2. Mask perturbation genes from response features.
3. Split by connected perturbation-gene component.
4. Add the public Kemmeren deletion-expression compendium
   [NCBI GEO GSE42528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42528)
   as the expression-matched yeast arm; GEO lists 2,633 *S. cerevisiae* samples.
   Retain SGA as an independent phenotype arm.
5. Only after the bridge passes, attach HIP/HOP for prospective compound ranking
   and measured turnaround-time validation.
