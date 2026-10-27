# Rebuilding the six manuscript figures

Download `YeastBridge-paper-figures-20261027.tar.gz` from the [extended release](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/paper-reproduction-20261027). It contains the complete copied R templates, frozen plotting tables, native panel PDFs, nine editable schematic scenes, Blender projects and molecular render assets. The wrapper verifies the archive and all 401 source files before writing output.

To assemble the original 51 native panels into six main figures:

```bash
python scripts/paper_figures.py --archive /path/to/YeastBridge-paper-figures-20261027.tar.gz \
  --output outputs/assembled-figures
```

To redraw the 42 data panels through the complete original templates, retaining the nine approved conceptual illustrations and molecular render assets:

```bash
python scripts/paper_figures.py --archive /path/to/YeastBridge-paper-figures-20261027.tar.gz \
  --output outputs/redrawn-figures --redraw --rscript /path/to/Rscript
```

The statistical plots use R; the evidence wheel retains its original Python source. This command consumes frozen plotting tables. Raw-data reconstruction and statistical tests are separate commands in [FULL_REPRODUCTION.md](FULL_REPRODUCTION.md). Molecular pictures and the conceptual illustrations are retained artwork; the command does not perform new Blender ray tracing or new molecular searches.

Both operations write a self-contained project under the chosen output directory. The six PDFs and previews are in `project/delivery/yeastbridge/figures/`. Original project and submission files are never modified. The PDF-only export avoids creating redundant 600-dpi PNG/TIFF copies; the copied export source also retains the original full-resolution mode.

The tested macOS environment is recorded in [figure_environment.json](../repro/paper-20261027/extended/figure_environment.json). Required Python packages are PyMuPDF, reportlab, Pillow, lxml, pandas, NumPy, matplotlib and pycirclize. R packages include ggplot2, cowplot, jsonlite, ggrepel, ggalluvial, ComplexHeatmap, circlize and scales; Arial and Helvetica were available during verification. Different font installations or rendering libraries can change pixels without changing plotted values.

Verification compared every panel raster at 72 dpi and the dimensions of all six pages. All 42 freshly drawn result panels matched the submitted panel pixels; all six figures assembled from the frozen native panels also matched. The wrapper writes a per-panel comparison and exits with a nonzero status if rendering differs. The [machine-readable result](../repro/paper-20261027/extended/figures.json) records the tested scope.
