# Verification of manuscript reproduction

The published model reproduces the cross-species screening chain, including fresh calibration and screen-wide permutations. The maintenance keeps the submitted checkpoint, candidate table and scientific settings fixed. The table below separates reconstruction from frozen inputs, fresh calculations, and observed differences.

| Component | What was actually run | Result and record |
| --- | --- | --- |
| Public raw-data access | Downloaded the GSE125162 count file; downloaded the first E-MTAB-2391 archive | Count SHA-256 matched; 216 CEL files checked from 1 of 18 archives. [Counts](public_counts_download.json), [CEL access](public_cel_download.json) |
| Raw barcode processing | Reparsed all 3,850 CEL files already present on the server, without cached tag measurements | Response archive byte-identical. This does not claim a fresh download of all 18 archives. [Record](raw_cel.json) |
| Training-input preparation | Rebuilt counts, protein matrix and cell roster | All three file hashes matched. [Record](training_input_rebuild.json) |
| Protein embeddings | Recomputed all 8,727 original sequence records with ESM-2 | All three arrays and all three indices byte-identical. [Record](embeddings.json) |
| Full screening entry | Regenerated profiles from the paper checkpoint; recalculated both 1,000-permutation calibrations and both 100,000-permutation maximum nulls | 5,216 complete-profile associations; same 25 intersecting associations, seven targets and 18 compounds. Candidate probabilities agree to 10⁻¹². [Whole entry](end_to_end.json), [calibration](calibration.json), [fresh nulls](fresh_screen.json) |
| Functional and backbone comparisons | Reran every internal/external functional query, B2/GO paired retrieval and backbone comparison on their frozen inputs | Exact row-level agreement and exact backbone JSON. [Record](core_comparisons.json) |
| Representation grid | Re-evaluated the four scFoundation/scYeast combinations using archived representations | All four cells exact. Foundation-model pretraining and feature extraction were not repeated. [Record](representation_grid.json) |
| A2/B2/C2/scYeast route evaluation | Retrained T2–T5 evaluation models on fixed representations and data splits | B2 T3/T5 detail tables agree within 10⁻⁶; B2 T2 differs slightly. [Record](route_comparisons.json) |
| Graph comparator | Reran the registered T2/T4 protocols on the released SGA edge table | All recorded metrics exact. [Record](graph_comparison.json) |
| Structural evidence from paper poses | Recomputed six contact/model analysis outputs and all 25 prioritization statuses | Exact scientific agreement; six prioritized associations. Also passed in a newly created chemistry environment. [Original runtime](structural_evidence.json), [new environment](structural_clean_environment.json), [dependency fix](chemistry_clean_runtime.json) |
| Fresh structural searches | Executed 282 candidate/native jobs from prepared coordinates, with no archived poses or old cache | 222 comparable Vina affinities exact; GNINA scores differ. All 21 native-qualification decisions match. Missing per-seed references remain listed. [Search comparison](redocking.json), [public input coverage](public_docking_plan.json) |
| Prioritization from fresh poses | Recalculated contacts and models using new candidate/native poses, retaining frozen background controls | Five prioritized associations; YB-IC-001 changes because one new CLIC1 pose lacks the required PHE83.B contact. The screening candidates are unchanged. [Record](structural_fresh_poses.json), [prepared input/protocol check](fresh_pose_input_verification.json) |
| Figures | Redrew 42 result panels using the complete original templates; retained nine conceptual illustrations and molecular artwork | All 51 panel rasters agree; six assembled figures agree. [Panels](figures.json), [assembly](figure_assembly.json), [visual inspection](figure_visual_review.json) |
| Paper training | Completed all six epochs of the historical recipe in an isolated output directory, in 163.4 minutes | 129 of 268 tensors exactly equal; gene-table mean absolute error 0.0482887 and mean row cosine 0.9708302. The original weights remain the paper default. [Full comparison](paper_training_comparison.json), [run](paper_training_run.json), [epoch metrics](paper_training_epochs.json), [one-step probes](training_recipe_equivalence.json) |

## Observed differences

The completed training replay did not regenerate the submitted model exactly: 139 of 268 tensors differ, with no missing or extra tensors. The 6,733 × 768 gene table has a maximum absolute difference of 0.5387554. Its final validation loss was 0.1315. The September training log and exact runtime were not recovered, so the source of the discrepancy is not established; repeated one-step runs of the original script also varied in the tested runtime. [All observed training packages](training_requirements.txt).

B2 T2 AUPRC was 0.5612326 in the recorded replay versus 0.5621393 in the reference; AUROC was 0.8197449 versus 0.8194967. The check reports these differences and returns a nonzero status. No threshold or reference table was changed to make it pass.

The archived poses reproduce all six structural prioritizations. In the new GNINA search, YB-IC-001 retains a median CNN affinity above the comparator, but seed 20261030 loses the PHE83.B contact. The unchanged requirement for both hotspot contacts in every seed therefore yields five prioritized associations. This is a computational replay of selected candidate/native searches; the complete background search panel was not regenerated.

These records support reconstruction of the submitted results from the frozen paper assets. They do not establish bitwise regeneration of every result from a newly trained model and newly searched molecular poses. Original model weights and results remain available alongside the replay evidence.

## Environments and checks

Screening and comparison checks also passed in a newly created Python 3.10 environment installed from the pinned requirements using cached wheels, with no inherited site packages. The full screen again produced the same 25 associations, and all core comparisons agreed exactly: [screen](core_clean_screen.json), [comparisons](core_clean_comparisons.json), [environment and installation fix](core_clean_environment.json), [all installed packages](core_clean_requirements.txt). The three active comparison configurations now read the release-pinned pair roster, removing their implicit dependency on an old result directory.

Training, chemistry and figure rendering use their separately documented environments. The chemistry installation exposed a missing `pkg_resources` dependency; the corrected environment pins `setuptools==80.10.2`. The [figure runtime](figure_environment.json) records the R, Python and font context.

The automated tests cover asset integrity and overwrite protection, assay-preserving permutations, family projection, non-estimable residuals, structural score/unit comparison and missing docking references. The reproduction commands and expected return codes are in [FULL_REPRODUCTION.md](../../../docs/FULL_REPRODUCTION.md).
