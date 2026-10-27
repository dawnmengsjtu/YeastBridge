# Manuscript data sources

The target sequence supplies the human query; reusable yeast perturbation measurements supply the experimental response collection. Exact runtime versions are pinned in [paper/assets.json](paper/assets.json).

| Source | Role in the manuscript |
| --- | --- |
| [Lee HIP/HOP, E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391) | Chemical-genetic responses reconstructed using published 20-mer strain barcodes and array tag-probe mappings |
| [GSE125162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125162) | Yeast single-cell training corpus |
| [Kemmeren GSE42528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42528), processed matrix [GSE42526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42526) | Deletion-expression reference and perturbation-response evaluation |
| [Norman GSE133344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344) | Human perturbation resource and internal functional queries; component accessions remain in the original provenance records |
| [Replogle K562, Figshare v1](https://doi.org/10.25452/figshare.plus.20029387.v1) | Human reference profiles; normalized bulk file 35773217 |
| [UniProt](https://www.uniprot.org/), [SGD](https://www.yeastgenome.org/) | Sequences and identifiers used for frozen ESM-2 embeddings |
| [Gene Ontology](https://geneontology.org/), OrthoDB, OMA and InParanoiDB | Functional coordinates and cross-species relationships |
| Drug Repurposing Hub, 18 August 2025 tables; ChEMBL 37; PubChem; BindingDB | Frozen known-target comparisons, activity records and chemical identities |
| Piotrowski et al. corrected supplementary data | Independent chemical-genetic measurements, with exact chemical-identity matching |
| Experimental receptor structures, binding studies, Vina and GNINA records | Structural prioritization and pocket-reference comparisons |

## Current processed response

`raw/paper/strain_response.npz` has SHA-256 `fef25d836b76d7654fd580e67fe316bb002268cc7fa94849b7cf754177c58471`. It contains 5,869 strains and 3,850 condition columns. Screening uses 5,865 represented, measured strains and 3,356 compound conditions corresponding to 3,250 compounds. UP/DOWN barcode contrasts are log2 relative abundance; the retained array key `z_score` does not imply a z statistic. The published strain file preserves HIP/HOP membership for within-assay permutation.

The prior ORF-probe response and competition matrix remain historical inputs and are not selected by the paper configuration. Frozen-pair availability yields 752 usable pairs: 510 train, 100 validation and 142 test. The 753-row pre-availability roster remains byte-identical in the inference archive.

## Availability and rights

The paper release includes the processed runtime inputs and the exact manuscript Supplementary Data 1. Raw public datasets can be obtained through their accessions; source files and original provenance keep their recorded dates. [Third-party terms](THIRD_PARTY_DATA.md) and [model terms](../models/training/vendor/scfoundation/MODEL_LICENSE) remain applicable. No new human or animal experiments were conducted for this project.
