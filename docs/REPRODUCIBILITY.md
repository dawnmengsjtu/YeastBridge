# Reproducing the manuscript

The current paper follows human target sequence → yeast strain profile → measured compound responses → screen-wide and family-residual testing → structural prioritization. Its assets are pinned independently of the historical competition release.

## Install and export

Use Python 3.9–3.12. Install `requirements.txt` for table export and numerical reconstruction; add `requirements-paper-model.txt` for frozen-model profile generation. The reference table is part of Git and needs no large download.

```bash
python -m pip install -r requirements.txt
python main.py --mode check --level reference
python predict.py --output outputs/published/results.csv
```

This exports 25 associations, seven targets and 18 compounds. `results.metadata.json` identifies this as a stored-result export. Changing the working directory does not change the chosen source.

## Runtime archives

Obtain archives from the [paper release](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/paper-ion-channel-20261027). The release checksum list verifies archives; `data/paper/assets.json` verifies every extracted file. Archives do not replace different existing assets. They contain regular, registered files only, with no links outside the checkout.

```bash
python main.py --mode install --archive /path/to/YeastBridge-paper-numerical-20261027.tar.gz
python main.py --mode install --archive /path/to/YeastBridge-paper-inference-20261027.tar.gz
python main.py --mode check --level all
```

| Archive | Contents | Purpose |
| --- | --- | --- |
| `YeastBridge-paper-inference-20261027.tar.gz` | Original B2 checkpoint and gene table; protein embeddings; 753-row frozen pair roster; yeast metadata; corrected response matrix; assay membership and compound identities | Rebuild the frozen model's human-to-yeast profiles |
| `YeastBridge-paper-numerical-20261027.tar.gz` | Complete/residual standardized rank inputs, observed/calibration arrays, both sets of 100,000 maximum-null values and their plans | Recompute observed scores and screen-wide decisions |
| `YeastBridge-paper-source-data-20261027.zip` | The manuscript's Supplementary Data 1, including figure source tables, analysis scripts/configurations, structural coordinates and numerical records | Trace manuscript figures and supplementary analyses to their inputs |

The inference archive is about 567 MB compressed; the numerical archive is about 197 MB. Allow at least 4 GB free disk for archives, verified installation and profile outputs. Installer staging also requires free temporary space. The numerical check and profile generation run on CPU. Historical full-null and structural jobs used separate GPU/chemistry environments recorded with their analyses.

## Numerical reconstruction

```bash
python main.py --mode reproduce --output outputs/numerical
```

This independently multiplies the rank-standardized target and condition matrices in float64, repeats condition maximization for both directions, and compares every observed correlation with the frozen scores (absolute tolerance 5×10⁻⁷). It then applies the add-one screen-wide estimator to the exact archived observed/calibration arrays and the saved maximum-null distributions. This preserves the original decision values despite CPU/GPU rounding.

The complete layer covers 2,197,000 directional target–compound hypotheses; the residual layer covers 617,500. The expected output is 5,216 complete-profile associations and 25 associations in the two-test intersection. Candidate identities, probabilities and correlations are compared with the reference table. `verification.json` records the actual errors and explicitly records that fresh permutations and training were not run.

## Rebuild profiles from the frozen model

```bash
python -m pip install -r requirements-paper-model.txt
python main.py --mode profiles --output outputs/profiles
```

The command uses `configs/paper_tasks.json` and `configs/paper_assets.json`. It verifies each input; fits the three frozen alignment operators on the same 510 training pairs; exports RRF profiles over 5,865 strains; selects the 338-entry target panel; and recalculates its common-mode correction. It checks the reference variance components and, when numerical assets are installed, compares the entire rank-standardized 338×5,865 matrix with the frozen screening input. The original checkpoint, raw-before-projection B2 input, member definitions and k=20 are retained.

Each run writes to a new directory. It does not overwrite frozen inputs or silently substitute stored candidate tables if computation fails.

## Scope of reproduction

| Operation | Current support |
| --- | --- |
| Frozen candidate export and reference checksums | Executable from Git alone |
| Every observed complete/residual correlation, screen-wide P and candidate intersection | Executable with the numerical archive |
| Human protein representation to corrected yeast profiles | Executable with the inference archive; rank comparison additionally uses numerical inputs |
| Retraining the exact paper checkpoint from initial weights | Historical code/evidence supplied; exact end-to-end retraining replay is not established by this maintenance |
| All raw-CEL preprocessing | Public accessions and scripts supplied; full raw downloads and platform-dependent preprocessing are separate |
| Regenerating both 100,000-permutation maximum nulls | Original GPU analysis code/settings supplied; current reconstruction uses their frozen distributions |
| All model comparisons, docking, structural refinement and figure rendering | Source records supplied; not all are integrated into a single portable rerun command |

Accordingly, a clone plus the release assets supports the specific computational reconstruction above. It does not establish a fresh end-to-end rerun of every experiment. Historical competition training and screening commands are described in their dated release, not used as substitutes for the manuscript.

## Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover export provenance and relocation, checksum failures, rejection of unregistered archive contents, finite Monte Carlo estimates and separate directional decisions. Run records from the maintenance are in `repro/paper-20261027/`.
