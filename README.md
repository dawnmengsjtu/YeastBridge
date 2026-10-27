# YeastBridge

**YeastBridge enables cross-species drug screening for human ion channels using yeast perturbation profiles**

YeastBridge connects human target sequences to measured yeast perturbation responses. A human target is represented as a profile across yeast strains, which queries an existing compound-response collection to prioritize compounds. This cross-species connection makes low-cost yeast measurements reusable across human targets without new human screening experiments.

The framework was selected through comparisons of GO and ESM-2 representations, Geneformer/scGPT/scFoundation backbones, and protein injection, ortholog mapping, training from scratch, scYeast and graph propagation. The final route combines direct sequence matching, supervised alignment and scFoundation protein injection using reciprocal rank fusion (`k=20`), followed by shared-component correction and compound matching.

[Current candidates](results.csv) · [Method](docs/CURRENT_METHOD.md) · [Reproduction guide](docs/REPRODUCIBILITY.md) · [Model](models/MODEL_CARD.md) · [Data](data/DATA_SOURCES.md)

## Manuscript results

The results follow three questions: model performance, model composition and compound screening.

- **Performance:** joint human–yeast functional retrieval achieves pooled MRR 0.4004, 9.1% above human references. Training-data experiments retain most of the full-data response-prediction performance with one-quarter of the training perturbations.
- **Composition:** sequence-based matching and scFoundation protein injection connect the two species; the comparison includes GO coordinates, other transcriptomic backbones and yeast baselines.
- **Screening:** 338 ion-channel-related entries and 3,250 compounds across 5,865 represented, measured yeast strains yield **25 candidate associations across seven targets and 18 compounds**. Structural analyses prioritize six associations. Candidates and functional retrieval are the reported endpoints; the structural scores are computational evidence.

## Quick start

```bash
git clone https://github.com/dawnmengsjtu/YeastBridge.git
cd YeastBridge
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python main.py --mode check
python predict.py --output outputs/published/results.csv
```

Python 3.9–3.12 is supported by the paper dependency set. The export writes the frozen 25-row table and a provenance record; it does not rerun screening. Paths resolve from the checkout, so commands can be invoked from another working directory.

## Reconstruct the analysis

Download the `inference` and/or `numerical` archives from the [paper release](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/paper-ion-channel-20261027). Each installed file is checked against [the immutable asset manifest](data/paper/assets.json).

```bash
python main.py --mode install --archive /path/to/YeastBridge-paper-numerical-20261027.tar.gz
python main.py --mode reproduce --output outputs/numerical

python -m pip install -r requirements-paper-model.txt
python main.py --mode install --archive /path/to/YeastBridge-paper-inference-20261027.tar.gz
python main.py --mode profiles --output outputs/profiles
```

`reproduce` independently recomputes all observed correlations for the complete and residual screens, then reconstructs the screen-wide probabilities and candidate intersection using the frozen permutation calibration and maxima. `profiles` regenerates the human-sequence-to-yeast profiles using the paper checkpoint, restricts common-mode correction to the 338-entry panel, and checks the resulting rank matrix when the numerical archive is installed. Both require a new output directory.

The [reproduction guide](docs/REPRODUCIBILITY.md) distinguishes these verified operations from retraining, regenerating the 100,000-permutation nulls, complete raw-CEL preprocessing and all structural analyses. The repository does not claim that a clone alone reruns every manuscript experiment.

## Versions and files

The manuscript model is pinned by SHA-256 `112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7`. The separately trained `submission-v2-20261003` checkpoint and 80-pair competition workflow remain historical releases; they are not the paper defaults. [Version map](docs/VERSIONS.md).

| Path | Role |
| --- | --- |
| `main.py`, `predict.py`, `run.sh` | Current manuscript entry points |
| `configs/paper_*.json` | Paper checkpoint, representations and frozen matching settings |
| `data/paper/`, `results.csv` | Current target roster, candidates, evidence and asset checksums |
| `scripts/` | Matching, screening, evaluation and reconstruction code |
| `models/` | Model provenance, training implementations and third-party licences |
| `docs/` | Current guides and explicitly dated method records |
| `raw/paper/` | Installed runtime data; supplied as release assets |
| `repro/`, older `results_a6000/`, `tasks/` | Dated historical calculations |

Validation: `python -m unittest discover -s tests -v`. Third-party conditions remain in [NOTICE.md](NOTICE.md) and the [data licence register](data/THIRD_PARTY_DATA.md). Model-weight licensing is retained in the repository. Citation metadata are in [CITATION.cff](CITATION.cff).
