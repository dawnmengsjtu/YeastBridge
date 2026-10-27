# Extended manuscript reproduction

The central executable chain is human target → transferred yeast profile → measured compound response → complete-profile and family-residual tests → candidate intersection. All commands below write to new directories and retain the published checkpoint, target panel, 5,865-strain universe, 510 fitting pairs, three RRF members and k=20.

Install the two core archives in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). Additional assets are in the [extended release](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/paper-reproduction-20261027):

```bash
python main.py --mode install --archive /path/to/YeastBridge-paper-preprocessing-20261027.tar.gz
python main.py --mode install --archive /path/to/YeastBridge-paper-comparisons-20261027.tar.gz
python main.py --mode check --level preprocessing
python main.py --mode check --level comparisons
```

The preprocessing archive contains the barcode-probe map, original SDRF and gene master. The comparison archive contains frozen benchmark inputs, representations, splits and reference results. Per-file size and SHA-256 checks are in `data/paper/extended_assets.json`. It includes the original A2/B2/C2 route-selection gene tables; their hashes are distinct from the later paper inference checkpoint and are checked before evaluation. Direct ESM-2 is not restored as an independent comparator.

The same release provides `YeastBridge-paper-code-20261027.tar.gz`, a compact executable source snapshot with a code-version record and file checksums. It omits large historical result directories; install the runtime archives above and the core archives after extraction. `YeastBridge-paper-verification-20261027.tar.gz` contains replay records, selected computed tables, fresh docking poses and training logs. It does not replace the submitted model with the newly trained weights.

## Raw measurements

The public [E-MTAB-2391 accession](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391) provides 3,850 CEL arrays in 18 archives. The downloader verifies every decompressed CEL against the exact file manifest, preserves different existing files, and writes a download report. Allow at least 20 GB for downloads, extracted CEL files and preprocessing outputs.

```bash
python scripts/download_paper_raw.py --dataset cel --output outputs/raw-cel
python scripts/recalc_fd_v1.py --build-response \
  --data-root outputs/raw-cel \
  --barcode-map raw/preprocessing/barcode_probe_map.tsv \
  --strain-universe raw/paper/published_strain_universe.tsv \
  --workers 8 --out outputs/rebuilt-response
```

A fresh parse of all 3,850 arrays reproduced the response archive byte for byte: SHA-256 `fef25d836b76d7654fd580e67fe316bb002268cc7fa94849b7cf754177c58471`. No cached tag measurements entered that run. The response contains 5,869 mapped strains; the model represents 5,865 of them.

GSE125162 counts can be downloaded and converted separately:

```bash
python scripts/download_paper_raw.py --dataset counts --output outputs/raw-counts
python scripts/prepare_training_data.py \
  --counts-tsv outputs/raw-counts/GSE125162_ALL-fastqTomat0-Counts.tsv \
  --output outputs/rebuilt-training-inputs
```

All three rebuilt training inputs—counts, protein matrix and cell roster—match their frozen SHA-256 values. This operation prepares inputs; the training procedure is in [PAPER_TRAINING.md](PAPER_TRAINING.md).

## Protein sequence representations

The original FASTA inputs and checkpoint digests are in `data/provenance/esm2/`. In the training environment, install `fair-esm==2.0.0` and `pyarrow==17.0.0`. Obtain `esm2_t33_650M_UR50D.pt` and its adjacent `esm2_t33_650M_UR50D-contact-regression.pt` from the [official ESM model release](https://github.com/facebookresearch/esm#available-models). Both files are hash-checked by the entry:

```bash
python scripts/paper_embeddings.py --model /path/to/esm2_t33_650M_UR50D.pt \
  --output outputs/rebuilt-embeddings --device cuda
```

This rebuilds all 8,727 original sequence records: 6,742 yeast proteins, 810 human alignment records and 1,175 records in the historical screening universe. The paper applies its fixed 338-entry target panel downstream. The fresh replay matched the SHA-256 of all three index files and all three float32 embedding arrays. The ESM representation is shared by the final fusion members; this entry does not introduce an independent ESM-2 comparator.

## Fresh screening

Use the paper-model environment and a CUDA-capable GPU. The tested configuration is Python 3.9, PyTorch 2.4.0+cu121, NumPy 1.26.4, pandas 2.2.2, SciPy 1.11.4, an RTX A6000 and TF32 disabled.

```bash
python scripts/run_paper_full.py \
  --response outputs/rebuilt-response/strain_response.npz \
  --output outputs/fresh-screen --device cuda:0
```

The command regenerates profiles, applies the frozen family projection, builds both rank matrices, calculates each 1,000-permutation calibration distribution with seed 42, then generates both 100,000-permutation maximum distributions with seed 20261024. Every permutation stays within HIP/HOP assay strata and repeats dose and direction maximization. It checks the complete 338-target and estimable 95-target residual universes before forming the intersection.

The stages can also be run separately:

```bash
python main.py --mode profiles --output outputs/profiles
python scripts/paper_calibration.py --profiles outputs/profiles/corrected_profiles \
  --response outputs/rebuilt-response/strain_response.npz --output outputs/calibration
python scripts/paper_permutations.py --calibration outputs/calibration --output outputs/permutations
```

The expected final result is 5,216 complete-profile associations and 25 intersecting associations across seven targets and 18 compounds. The checker retains full-precision newly computed arrays, allows 5×10⁻⁷ absolute correlation roundoff, and checks candidate probabilities to 10⁻¹². Frozen TSV serialization and CPU/GPU reductions can change low-order correlation digits; the reported errors remain in the verification file.

## Model comparisons

Install `requirements-paper-comparisons.txt`. The core entry checks every query/paired-retrieval row and the complete backbone result, including internal/external functional retrieval, B2 and GO alignment, and Geneformer/scGPT/scFoundation selection.

```bash
python -m pip install -r requirements-paper-comparisons.txt
python scripts/paper_comparisons.py --mode core --output outputs/comparisons
python scripts/paper_comparisons.py --mode representation-grid --output outputs/representation-grid
python scripts/paper_comparisons.py --mode routes --output outputs/routes
```

The representation grid uses archived scFoundation/scYeast features and a verified K562 observation index; it evaluates their combinations without repeating foundation-model pretraining. The route entry retrains T2–T5 evaluation models on frozen A2/B2/C2/scYeast representations and splits. B2 T3 and T5 are checked against the registered B2 detail tables. Graph propagation has its own registered T2/T4 protocols and separate entry:

```bash
python main.py --mode install --archive /path/to/YeastBridge-paper-graph-20261027.tar.gz
python scripts/paper_graph_comparison.py --output outputs/graph-comparison
```

The route checker reports numerical differences rather than overwriting reference values. The current B2 T2 logistic-regression replay differs slightly from its reference; thread-count probes also varied, but the exact cause of the reference discrepancy has not been established. The recorded differences and runtime remain available. T3/T5 detailed checks use a 10⁻⁶ tolerance. A nonzero return from the route comparison indicates a recorded mismatch and must be inspected.

## Structural evidence

Download the unchanged `YeastBridge-paper-source-data-20261027.zip` from the core paper release. Create the chemistry environment with [environments-paper-structural.yml](../environments-paper-structural.yml):

```bash
conda env create -f environments-paper-structural.yml
conda activate yeastbridge-paper-structural
```

A newly created environment reproduced all six structural-analysis outputs and all 25 prioritization statuses. The initial installation exposed a missing `pkg_resources` dependency; the tested `setuptools==80.10.2` pin is now included in the YAML. The [runtime record](../repro/paper-20261027/extended/chemistry_clean_runtime.json) records the installation and fix. Install [GNINA 1.3.2](https://github.com/gnina/gnina/releases/tag/v1.3.2) using the default `gnina.1.3.2` release asset for fresh searches. The entry checks its SHA-256 against the recorded binary. In the chemistry environment, run:

```bash
python scripts/paper_structural_evidence.py \
  --archive /path/to/YeastBridge-paper-source-data-20261027.zip \
  --output outputs/structural-evidence
```

This verifies the source ZIP, creates an isolated working copy, recalculates residue contacts, fits the contact models, and reapplies the manuscript's descriptive prioritization. It compares the scientific content of six analysis outputs and all 25 candidate statuses. The expected prioritized IDs are YB-IC-001, 002, 021, 022, 023 and 025. The descriptor is separate from the stricter binding-control classification stored in `pair_verdicts.tsv`; neither represents experimental validation of candidates.

For fresh candidate/native searches from prepared coordinates, extract `structural_current/` from the source ZIP and run:

```bash
python scripts/paper_redocking.py --source /path/to/structural_current \
  --gnina /path/to/gnina --gnina-devices 0 --workers 8 \
  --output outputs/redocking
```

The entry retains receptor, ligand, pocket, search-depth and seed settings. It copies no previous poses or docking cache. A custom runtime launcher can be supplied as a JSON array with `--gnina-launcher-json`. GNINA 1.3.2, Vina and the recorded chemistry dependencies are required; the [binary digest and observed launcher](../repro/paper-20261027/extended/redocking_runtime.json) identify the executable and the compatibility loader used on the older server; see the [observed chemistry environment](../repro/paper-20261027/extended/structural_environment.json). Complete background-panel searches and every historical structural variant are separate from this candidate/native audit.

`--plan-only` prepares the inputs and checks their coverage without running searches. The public source ZIP produced the same 282 job identities as the completed replay. All 282 jobs completed; 222 comparable Vina scores matched exactly, while 12 comparable GNINA scores differed. The archive supplies per-seed native scores for 15 jobs and aggregate qualification records for the remaining native cases; all 21 available native-qualification decisions matched. The full comparison keeps missing per-seed references explicit.

To propagate the newly searched poses through contact analysis and prioritization:

```bash
python scripts/paper_structural_evidence.py \
  --archive /path/to/YeastBridge-paper-source-data-20261027.zip \
  --redocking outputs/redocking --output outputs/structural-fresh-poses
```

This audit retains frozen background-control searches and replaces the candidate/native jobs. It verifies the required prepared inputs, saves copies of its replay inputs and stops if native qualification changes before reusing pair-verdict metadata. The recorded replay prioritized five associations: YB-IC-002, 021, 022, 023 and 025. YB-IC-001 lost the PHE83.B contact in one of its three new CLIC1 poses and therefore did not satisfy the unchanged all-seed contact rule. Its median CNN affinity remained above the comparator. The checker returns a nonzero status for this changed prioritization; the published six-association result remains the result of the archived poses.

## Figures

The [figure reproduction guide](FIGURE_REPRODUCTION.md) covers assembly of all six figures and template-based redrawing of 42 result panels. The complete plotting project is provided as a separate verified archive.

## Verification records

The [verification ledger](../repro/paper-20261027/extended/README.md) and its machine-readable records distinguish exact equality, numerical tolerance, fresh computations and reused inputs. Original paper results are not replaced by replay outputs. The core tests cover installer integrity and Monte Carlo decision logic; extended tests cover assay-preserving permutations, family projection and zero-residual handling.
