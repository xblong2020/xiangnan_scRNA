suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript export_seurat_counts_to_mtx.R <input.rds> <output_dir>")
}

input_rds <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
assay <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay, layer = "counts")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

matrix_path <- file.path(output_dir, "counts.mtx")
barcodes_path <- file.path(output_dir, "barcodes.tsv")
features_path <- file.path(output_dir, "features.tsv")
obs_path <- file.path(output_dir, "obs.tsv")

writeMM(counts, matrix_path)
writeLines(colnames(counts), barcodes_path, useBytes = TRUE)
writeLines(rownames(counts), features_path, useBytes = TRUE)

obs <- obj@meta.data[colnames(counts), , drop = FALSE]
write.table(obs, obs_path, sep = "\t", quote = FALSE, col.names = NA)

cat("WROTE", matrix_path, "\n")
cat("WROTE", barcodes_path, "\n")
cat("WROTE", features_path, "\n")
cat("WROTE", obs_path, "\n")
cat("assay", assay, "features", nrow(counts), "cells", ncol(counts), "nnz", length(counts@x), "\n")
