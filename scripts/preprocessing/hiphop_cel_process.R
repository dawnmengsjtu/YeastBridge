#!/usr/bin/env Rscript
# HIP/HOP E-MTAB-2391 strain-barcode expression from raw Affymetrix CEL files.
#
# Probe-set structure comes from the platform TagProbes.tab coordinates
# (probeset, x, y) joined by probe sequence - never from guessed index
# columns or binary CDF parsing. Cell order follows the standard Affymetrix
# layout: linear_index = x + y * n_cols (0-based).
#
# Method (declared exactly): affyio cell intensities -> log2 -> quantile
# normalisation across arrays -> median-polish summarisation per probe set.
# This is RMA without the convolution background model; no stronger method
# claim is made anywhere downstream.
#
# All paths are arguments; nothing is hardcoded.

suppressMessages({
  library(affyio)
  library(preprocessCore)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: hiphop_cel_process.R <cel_list.tsv> <probeset_map_xy.tsv> <out_expr.gz> <out_summary.txt>")
}
cel_list_path <- args[1]
map_path <- args[2]
out_expr <- args[3]
out_summary <- args[4]

cel_table <- read.delim(cel_list_path, header = FALSE,
                        col.names = c("array", "path"),
                        stringsAsFactors = FALSE)
cat(sprintf("arrays: %d\n", nrow(cel_table)))

header <- read.celfile.header(cel_table$path[1])
dims <- header$`CEL dimensions`
n_cols <- as.integer(dims[["Cols"]])
n_cells_total <- n_cols * as.integer(dims[["Rows"]])
cat(sprintf("chip: %s  %d cols\n", header$cdfName, n_cols))

map <- read.delim(map_path, header = TRUE,
                  col.names = c("probeset", "x", "y", "probe"),
                  stringsAsFactors = FALSE)
indices <- map$x + map$y * n_cols + 1L
stopifnot(all(indices >= 1 & indices <= n_cells_total))
probe_sets <- split(indices, map$probeset)
cat(sprintf("probe sets: %d  cells: %d\n",
            length(probe_sets), length(indices)))

read_cells <- function(path) {
  ab <- read.celfile(path)
  ab$INTENSITY$MEAN
}

first <- read_cells(cel_table$path[1])
stopifnot(length(first) == n_cells_total)
intensities <- matrix(NA_real_, nrow = n_cells_total, ncol = nrow(cel_table))
colnames(intensities) <- cel_table$array
intensities[, 1] <- first
for (i in seq_len(nrow(cel_table))[-1]) {
  intensities[, i] <- read_cells(cel_table$path[i])
  if (i %% 100 == 0) cat(sprintf("  read %d/%d\n", i, nrow(cel_table)))
}

expr <- log2(intensities)

# Quantile normalisation in pure R (same algorithm as preprocessCore,
# which cannot spawn its worker threads in every containerised runtime):
# target distribution = column-mean of per-column sorted values.
sorted <- apply(expr, 2, sort)
target <- rowMeans(sorted)
expr_norm <- vapply(seq_len(ncol(expr)), function(j) {
  target[rank(expr[, j], ties.method = "average")]
}, numeric(n_cells_total))
rm(sorted, intensities)

summarise_one <- function(idx) {
  block <- expr_norm[idx, , drop = FALSE]
  fit <- medpolish(block, trace.iter = FALSE, maxiter = 50)
  fit$overall + fit$col
}

set_ids <- names(probe_sets)
result <- t(vapply(set_ids, function(ps) {
  summarise_one(as.integer(probe_sets[[ps]]))
}, numeric(ncol(expr))))
colnames(result) <- cel_table$array

con <- gzfile(out_expr, "w")
write.table(result, con, sep = "\t", quote = FALSE)
close(con)

lines <- c(
  sprintf("arrays_processed\t%d", nrow(cel_table)),
  sprintf("cells_per_array\t%d", n_cells_total),
  sprintf("chip_type\t%s", header$cdfName),
  sprintf("probe_sets\t%d", length(set_ids)),
  sprintf("mapped_cells\t%d", length(indices)),
  "method\tlog2+pure_R_quantile_normalisation+median_polish_no_bg_model"
)
writeLines(lines, out_summary)
cat("done\n")
