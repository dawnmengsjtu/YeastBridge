# Third-party sources and distribution notes

- ESM: https://github.com/facebookresearch/esm and its MIT LICENSE. The project reads precomputed ESM-2 embeddings in its full analysis mode.
- scFoundation: https://github.com/biomap-research/scFoundation. Source code is Apache-2.0; model weights have a separate non-commercial research MODEL_LICENSE. B2 was trained using this backbone. Preserve applicable upstream notices when distributing a checkpoint; this repository does not relabel it as unrestricted original weights.
- Historical scGPT, Geneformer and scYeast comparisons are identified in `docs/model_selection/`. They are not bundled third-party source/model distributions in the review package, and their old training environments have not been reconstructed here.
- Dataset accessions, retained artifacts and unresolved license records are listed in `data/DATA_SOURCES.md`. Public accessibility alone is not an assertion of redistribution permission.

This maintenance change does not assign a new project-wide license or replace upstream license terms. No external model weight is included in the review bundle. Numerical results remain subject to the statistical limitations described in the model card and protocols.
