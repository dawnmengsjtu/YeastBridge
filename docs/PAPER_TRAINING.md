# Replaying the paper training recipe

The paper inference checkpoint is fixed at SHA-256 `112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7`. Its recovered training recipe is `configs/b2_training.json`, implementation `historical-v1`. `train.py` defaults to that recipe. The later competition model uses the explicitly separate `configs/b2_training_v2.json`.

Use the training environment separately from the screening environment:

```bash
python3.10 -m venv .venv-training
source .venv-training/bin/activate
python -m pip install -r requirements-training.lock.txt
python scripts/install_training_assets.py --checkpoint /path/to/models.ckpt
python train.py --config configs/b2_training.json --device cuda --output outputs/paper-training
```

The initial scFoundation checkpoint comes from its [official model page](https://github.com/biomap-research/scFoundation/blob/main/model/README.md); its expected SHA-256 is `9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9`. The installer checks the initial weights and data. Upstream [MODEL_LICENSE](../models/training/vendor/scfoundation/MODEL_LICENSE) conditions apply. The [raw preprocessing guide](FULL_REPRODUCTION.md) also reconstructs all training inputs from the source counts and frozen protein embeddings.

The full run uses 38,225 cells, 6,733 genes, six epochs, batch 32, seed 42, learning rate 10⁻⁴, a 1,200-gene sequence cap and bf16. It preserves the original historical masking and validation behavior. The later v2 corrections are not silently applied to a paper replay. Each run saves the exact configuration, split, observed packages, log and output hashes.

A successful six-epoch run establishes execution of the recovered recipe. Exact regeneration of the September paper checkpoint must additionally be demonstrated by comparing learned tensors and the projected gene table, not merely by the existence of a new checkpoint or a similar loss curve. The contemporary September training log was not recovered; earlier August logs are retained as earlier records. The maintenance comparison reports the actual agreement obtained in the recorded runtime.

An identical 32-cell, one-step probe also varied when the original script itself was repeated: maximum gene-table absolute difference 7.31×10⁻⁵. The original-versus-portable probe differed by 7.99×10⁻⁵. These observations are recorded in [training_recipe_equivalence.json](../repro/paper-20261027/extended/training_recipe_equivalence.json); they establish numerical variation in the tested runtime, without identifying its exact cause or proving full-training equivalence.

After training, compare the learned arrays explicitly:

```bash
python scripts/compare_paper_training.py --run outputs/paper-training --output outputs/training-comparison
```

The checker compares every learned tensor, reports per-tensor absolute errors, and checks the exported gene table. It returns a nonzero status when tensors differ; matching architecture or a successful run is not treated as numerical equality.

The recorded six-epoch replay completed in 163.4 minutes with final validation loss 0.1315. Of 268 tensors, 129 were exactly equal and 139 differed, with no missing or extra tensors. The 6,733 × 768 gene table had mean absolute error 0.0482887, maximum absolute error 0.5387554 and mean row cosine 0.9708302 against the submitted model. The [full comparison](../repro/paper-20261027/extended/paper_training_comparison.json), [run record](../repro/paper-20261027/extended/paper_training_run.json) and [epoch metrics](../repro/paper-20261027/extended/paper_training_epochs.json) preserve this result. The exact cause of the discrepancy has not been established.

Use published weights for reconstruction of submitted results. A newly trained model is a separate result and does not automatically become the paper default.
