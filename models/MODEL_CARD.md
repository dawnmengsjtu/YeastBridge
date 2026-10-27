# Manuscript model

YeastBridge converts a human target sequence into a yeast-strain query and uses measured compound responses for human ion-channel compound prioritization. The cross-species interface and reuse of low-cost yeast measurements are the primary contributions.

## Frozen components

| Component | Paper configuration |
| --- | --- |
| Protein representation | ESM-2 `esm2_t33_650M_UR50D`, layer 33 residue-mean pooling, 1,280 dimensions |
| Protein injection | scFoundation-derived yeast model, raw ESM-2 → affine 1,280-to-768 projection → normalized comparison |
| Paired mapping | 510 training pairs; 3 removed PCs; ridge=1, anchored ridge=1, PLS=8 |
| Fusion | `direct_pc1`, `aligned_c_ensemble`, `b2_verbatim`; RRF k=20 |
| Screening panel | 338 ion-channel-related entries, 5,865 measured represented strains, 3,250 compounds |
| Candidate rule | Complete-profile and family-residual screen-wide P≤0.05 for the same target, compound and direction |

All protein representations can be rebuilt from the frozen FASTA inputs with `scripts/paper_embeddings.py`; the sequence and checkpoint hashes are checked before inference. See the [extended reproduction guide](../docs/FULL_REPRODUCTION.md).

The original paper checkpoint has SHA-256 `112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7`; its yeast gene table has SHA-256 `d5cd2912f7f798f236a41a31e2d23509402480415f74129e66a2c865dbae9c24`. The inference installer verifies both. Configuration: [paper_assets.json](../configs/paper_assets.json).

Geneformer, scGPT and scFoundation were compared before selecting the backbone; ortholog mapping, protein injection, training from scratch, scYeast and graph propagation informed yeast-transfer selection. The direct sequence member remains inside the final fusion, not an additional independent method-selection arm.

The separately trained competition checkpoint (`281d363f…09c27ad`) and its six-epoch logs remain in the historical release. Those logs are not substituted for provenance of the paper checkpoint. The current profile command uses the frozen paper weights; it does not claim an exact retraining replay.

Outputs are computational target–compound associations. Functional retrieval, yeast response prediction and structural prioritization retain their respective evaluation meanings.

## Third-party terms

ESM-2 code and pretrained model terms remain with its upstream release. scFoundation code is Apache-2.0; its model weights and derived B2 model retain the separate [non-commercial research model licence](training/vendor/scfoundation/MODEL_LICENSE). See [NOTICE.md](../NOTICE.md).
