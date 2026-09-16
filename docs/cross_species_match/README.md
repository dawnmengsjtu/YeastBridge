# Cross-species KO / perturbation matching benchmark

This module asks one narrow, falsifiable question:

> On a locked set of unseen human or mouse perturbations, does adding yeast
> deletion/HIP-HOP information improve functional matching over either source
> alone?

It does **not** preselect GPCRs, ion channels, diseases or candidate targets. The
full input panels are projected into an explicitly supplied conserved-function
space. The benchmark is intentionally small-model and offline so that code can be
submitted and audited before large datasets or foundation models are ready.

## One-command synthetic check

From the repository root:

```bash
python scripts/cross_species_match.py \
  --config configs/cross_species_match.synthetic.json
pytest -q tests/test_cross_species_match.py
```

The bundled synthetic example intentionally ends in
`NO_GO_JOINT_ROUTE_SINGLE`: the joint arm ties the yeast arm, so it is not legal
to claim that integration helped. This demonstrates the guardrail rather than a
manufactured win.

## Inputs

The JSON config points to five tab-separated files.

### Signed matrices

`higher_matrix.tsv` and `yeast_matrix.tsv` are wide matrices. The first column is
`signature_id`; every remaining column is a species-native gene. Values are signed
effects such as z-scores or log fold changes. Human and mouse rows may coexist in
the higher-organism matrix; genes irrelevant to a row are zero.

### Signature metadata

Both metadata files require:

| column | meaning |
|---|---|
| `signature_id` | globally unique row ID |
| `species` | e.g. `human`, `mouse`, `yeast` |
| `function_label` | benchmark truth such as a phenotype or conserved program; not a target filter |
| `perturbation_genes` | one or more KO/perturbed genes, separated by `+`, `;`, `,` or `|` |
| `split_group` | indivisible biological lineage for splitting |
| `split` | locked `train`/`test`, or blank for deterministic stratified assignment |
| `assay_type` | optional, e.g. `CRISPR_KO`, `deletion`, `HIPHOP` |

All replicates, doses, cell lines, guide combinations and orthologous versions of
the same perturbation must use the same `split_group`. The program additionally
rejects an exact species-specific gene combination or exact signature duplicated
across train and test.

### Explicit projection table

`projection.tsv` requires:

| column | meaning |
|---|---|
| `species` | namespace of `gene` |
| `gene` | source matrix gene |
| `module_id` | common conserved functional module |
| `evidence` | `ortholog`, `go`, or `conserved_module` |
| `weight` | positive mapping confidence |

Evidence weights are fixed in the JSON before test evaluation. Each gene's
outgoing weights are normalised, signed effects are aggregated per module, and
the resulting vector is L2-normalised. Mapping coverage is reported for every
signature; low-coverage rows fail rather than silently becoming zero vectors.

## Locked three-arm comparison

Every arm predicts the same held-out higher-organism signatures:

1. `human_only`: nearest functional centroid learned only from the configured
   human/mouse train rows (the output records the exact species scope).
2. `yeast_only`: nearest functional centroid learned only from yeast train rows.
3. `joint`: a predeclared fixed-weight blend of the two score matrices.

Accuracy, macro-F1, mean reciprocal rank and top-k accuracy are emitted. Inference
is paired at `split_group`, using a paired bootstrap confidence interval and an
exact or Monte-Carlo sign-flip permutation test.

The prediction class universe is fixed without reading test labels: either list
`protocol.benchmark_function_labels` in config or let the program take the full
intersection of labels represented in both locked training sources. Train-only
distractor classes are retained.

`GO_JOINT` is allowed only when joint exceeds **both** single-source arms by the
predeclared minimum gain, both bootstrap lower bounds clear that gain, both
one-sided paired p-values clear alpha, and the minimum independent-group count is
met. Otherwise the run succeeds with `NO_GO_JOINT_ROUTE_SINGLE` and routes to the
stronger single arm. A tie break is fixed in config; the example chooses yeast for
its practical screening throughput, not because speed proves biological validity.

## Outputs

The versioned output directory contains:

- `metrics.json`: three-arm metrics, paired tests, decision, input hashes and audits.
- `predictions.tsv`: per-query paired predictions and ranks.
- `matches.tsv`: explicit top yeast deletion/HIP-HOP matches and module contributions.
- `projection_coverage.tsv`: per-signature mapping coverage.
- `split_lock.tsv`: immutable split assignment; changed reruns are rejected.
- `protocol_lock.json`: immutable code/config/input fingerprints and statistical gate.
- `RUN_SUMMARY.md`: short human-readable result.

Use a new output directory for a new dataset or protocol version. Do not overwrite
a lock after seeing test results.

## Claims this benchmark does not establish

Similarity is not drug efficacy. A passing computational bridge still needs a
prospective wet-lab panel with measured turnaround time, cost, assay success rate,
directional concordance and failure cases. Until then, say “prioritises a yeast
assay route” rather than “validates a human drug in yeast.”
