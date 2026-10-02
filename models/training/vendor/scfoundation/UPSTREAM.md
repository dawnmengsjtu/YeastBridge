# scFoundation source provenance

These model files were recovered from the original project environment on 2026-10-03 (`yeastbridge/src/external/scfoundation/model`). Their upstream Git commit was not recorded in that installation. The upstream project is https://github.com/biomap-research/scFoundation.

The source is covered by the included Apache-2.0 `LICENSE`; the separate `MODEL_LICENSE` applies to scFoundation model weights and their use. Both notices were retrieved from the upstream repository on 2026-10-03. B2 remains subject to the upstream weight terms, including the non-commercial research restriction.

Local changes to `load.py`: `load_model_frommmf` accepts an explicit device instead of always calling `.cuda()`; its trusted, checksum-verified checkpoint load explicitly uses `weights_only=False` for PyTorch 2.6 compatibility. The other vendored model files retain the recovered implementation. The project trainer loads only assets whose SHA-256 matches `configs/b2_training.json`.
