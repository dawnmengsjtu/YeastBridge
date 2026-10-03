# Third-party sources and distribution notes

- ESM-2: https://github.com/facebookresearch/esm, MIT license. The final pipeline reads frozen model outputs; the model name, extraction parameters, indices, recovered sequence inputs and checkpoint digest are documented in the model card and `data/provenance/esm2/`.
- scFoundation: https://github.com/biomap-research/scFoundation. Source code is Apache-2.0; model weights use the separate non-commercial research `MODEL_LICENSE`. Both notices accompany the vendored implementation. The B2 checkpoint is derived from this backbone and retains applicable upstream restrictions.
- UniProt/SGD protein and identifier data, orthology, GO annotations and expression datasets retain the terms and attribution listed in `data/THIRD_PARTY_DATA.md`. Exact content versions, use, source and acquisition evidence are recorded in `data/datasets.json`.
- scGPT, Geneformer and scYeast belong to historical model comparisons. They are not required by the current candidate-screening pipeline and are not redistributed as additional pretrained weights in the submission package.

This repository does not assign a new project-wide license or replace third-party terms. The review bundle contains no external model checkpoint; the full bundle contains the submitted B2 checkpoint and applicable model notice. Numerical results retain the scope and limitations stated in the model card.
