# Version map

| Version | Model / response | Result | Use |
| --- | --- | --- | --- |
| Paper ion-channel release, 2026-10-27 | B2 checkpoint `112b39af…5759f2be7`; corrected barcode response `fef25d83…77c58471` | 338 entries, 5,865 strains; 25 pairs / 7 targets / 18 compounds | Current manuscript and default entry points |
| Competition `submission-v2-20261003` | Separately trained checkpoint `281d363f…09c27ad`; earlier response matrix | 80 competition candidates | Historical release |
| September/early October mixed-target analyses | Older candidate panels and original coverage | Date-specific outputs | Historical methods and diagnostics |

File dates do not select a model. The checkpoint, processed-response checksum, target scope and permutation universe jointly identify a result. The paper configuration explicitly selects the manuscript model even though the competition checkpoint was trained later.

The active numerical sources are `data/paper/assets.json`, `configs/paper_tasks.json`, `configs/paper_assets.json` and the paper release archives. Historical records retain their original numbers and names; those numbers do not supersede the current manuscript.

The 2026-10-27 maintenance changes packaging, entry points, documentation and validation. It does not retrain weights, alter selection thresholds or replace scientific results.

`train.py`, `scripts/install_assets.py`, `scripts/check_assets.py`, `scripts/verify_submission.py`, `scripts/build_submission.py`, `models/submission.json`, `configs/b2_training_v2.json` and the `requirements-*-training/full` recipes belong to the preserved competition workflow. Current paper operations use `main.py`, `configs/paper_*.json`, `requirements.txt` and `requirements-paper-model.txt`. The original training sources are kept intact so recorded training-source hashes remain meaningful.
