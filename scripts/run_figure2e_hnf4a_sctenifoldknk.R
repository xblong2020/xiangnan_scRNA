#!/usr/bin/env Rscript

## HNF4A-only scTenifoldKnk virtual knockout in identity-high normal hepatocytes.
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
root <- normalizePath(if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
                      mustWork = FALSE)
if (!dir.exists(file.path(root, "scripts"))) root <- normalizePath(getwd(), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

input_dir <- get_arg("--input-dir",
  file.path(root, "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference"))
out_dir <- get_arg("--out-dir", file.path(root, "metadata/driver/figure2e_hnf4a_sctenifoldknk"))
target_tf <- get_arg("--target-tf", "HNF4A")
seed <- as.integer(get_arg("--seed", "11"))
nc_nnet <- as.integer(get_arg("--nc-nnet", "1"))
nc_ncells <- as.integer(get_arg("--nc-ncells", "100"))
nc_ncomp <- as.integer(get_arg("--nc-ncomp", "3"))
ma_ndim <- as.integer(get_arg("--ma-ndim", "2"))
ncores <- as.integer(get_arg("--ncores", "1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("R package scTenifoldKnk is required")
}
matrix_path <- file.path(input_dir, "figure2e_hnf4a_normal_reference_counts_genes_x_cells.mtx")
genes_path <- file.path(input_dir, "figure2e_hnf4a_normal_reference_genes.tsv")
counts <- as(Matrix::readMM(matrix_path), "dgCMatrix")
genes <- read.delim(genes_path, stringsAsFactors = FALSE)[[1]]
rownames(counts) <- make.unique(as.character(genes))
if (!target_tf %in% rownames(counts)) stop(target_tf, " is absent from the matrix")

set.seed(seed)
started <- Sys.time()
result <- scTenifoldKnk::scTenifoldKnk(
  countMatrix = counts, gKO = target_tf, qc = FALSE, qc_minCells = 3,
  nc_nNet = nc_nnet, nc_nCells = nc_ncells, nc_nComp = nc_ncomp,
  ma_nDim = ma_ndim, nCores = ncores
)
candidates <- c("diffRegulation", "diff_regulation", "diffRegulationResult")
perturb <- NULL
for (name in candidates) {
  if (!is.null(result[[name]]) && is.data.frame(result[[name]])) {
    perturb <- as.data.frame(result[[name]])
    break
  }
}
if (is.null(perturb)) stop("scTenifoldKnk result lacks diffRegulation table")
if (!"gene" %in% names(perturb)) perturb$gene <- rownames(perturb)
if (!"tf" %in% names(perturb)) perturb$tf <- target_tf
perturb$subset <- "normal_reference"
preferred <- c("tf", "gene", "distance", "Distance", "p.adj", "p_adj",
               "pvalue", "p.value", "FDR", "Z", "FC", "subset")
perturb <- perturb[, unique(c(intersect(preferred, names(perturb)),
                              setdiff(names(perturb), preferred))), drop = FALSE]
if ("Distance" %in% names(perturb) && !"distance" %in% names(perturb)) {
  names(perturb)[names(perturb) == "Distance"] <- "distance"
}
if ("p_adj" %in% names(perturb) && !"p.adj" %in% names(perturb)) {
  names(perturb)[names(perturb) == "p_adj"] <- "p.adj"
}
if ("FDR" %in% names(perturb) && !"p.adj" %in% names(perturb)) {
  names(perturb)[names(perturb) == "FDR"] <- "p.adj"
}
if ("pvalue" %in% names(perturb) && !"p.value" %in% names(perturb)) {
  names(perturb)[names(perturb) == "pvalue"] <- "p.value"
}

tsv_path <- file.path(out_dir, "figure2e_hnf4a_normal_reference_perturbation_genes.tsv")
rds_path <- file.path(out_dir, "figure2e_hnf4a_normal_reference_sctenifoldknk.rds")
report_path <- file.path(out_dir, "figure2e_hnf4a_normal_reference_sctenifoldknk_report.json")
write.table(perturb, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(result, rds_path)
report <- list(
  module = "Figure 2E HNF4A scTenifoldKnk", target_tf = target_tf,
  method = "scTenifoldKnk virtual knockout", subset = "normal_reference",
  n_genes = nrow(counts), n_cells = ncol(counts),
  parameters = list(qc = FALSE, qc_minCells = 3, nc_nNet = nc_nnet,
    nc_nCells = nc_ncells, nc_nComp = nc_ncomp, ma_nDim = ma_ndim,
    nCores = ncores, seed = seed),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  n_significant_excluding_target = sum(perturb$gene != target_tf &
    is.finite(as.numeric(perturb$p.adj)) & as.numeric(perturb$p.adj) < 0.05),
  outputs = list(perturbation_genes = normalizePath(tsv_path, winslash = "/", mustWork = FALSE),
                 result_rds = normalizePath(rds_path, winslash = "/", mustWork = FALSE),
                 report = normalizePath(report_path, winslash = "/", mustWork = FALSE)),
  caveat = "Network perturbation evidence from a virtual knockout; parameters match the existing Module 7 runs."
)
write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE)
message("HNF4A normal-reference scTenifoldKnk completed: ", tsv_path)
