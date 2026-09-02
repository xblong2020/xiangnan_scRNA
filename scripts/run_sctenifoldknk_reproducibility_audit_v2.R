#!/usr/bin/env Rscript

## Versioned three-axis scTenifoldKnk rerun for the Figure 2/3 reproducibility audit.

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  equal_hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(equal_hit)) return(sub(paste0("^", flag, "="), "", equal_hit[1]))
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

as_bool <- function(value) {
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

input_dir <- get_arg("--input-dir", stop("--input-dir is required"))
matrix_file <- get_arg("--matrix-file", "sctenifoldknk_counts_genes_x_cells.mtx")
genes_file <- get_arg("--genes-file", "sctenifoldknk_genes.tsv")
output_dir <- get_arg("--output-dir", stop("--output-dir is required"))
metadata_dir <- get_arg("--metadata-dir", stop("--metadata-dir is required"))
target_tf <- get_arg("--target-tf", stop("--target-tf is required"))
axis <- get_arg("--axis", "")
subset <- get_arg("--subset", "")
seed <- as.integer(get_arg("--seed", stop("--seed is required")))
nc_nnet <- as.integer(get_arg("--nc-nnet", "10"))
nc_ncells_requested <- as.integer(get_arg("--nc-ncells", "500"))
nc_ncomp <- as.integer(get_arg("--nc-ncomp", "3"))
ma_ndim <- as.integer(get_arg("--ma-ndim", "2"))
ncores <- as.integer(get_arg("--ncores", "8"))
qc <- as_bool(get_arg("--qc", "false"))
qc_min_cells <- as.integer(get_arg("--qc-min-cells", "3"))

path_for_report <- function(path) gsub("\\\\", "/", as.character(path))

if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("R package scTenifoldKnk is required")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
# Keep the ASCII subst path intact; normalizePath would expand it back to the
# non-ASCII OneDrive path and break R 4.5 file access on Windows.
matrix_path <- file.path(input_dir, matrix_file)
genes_path <- file.path(input_dir, genes_file)
if (!file.exists(matrix_path) || !file.exists(genes_path)) stop("Input matrix or gene index is missing")
counts <- as(Matrix::readMM(matrix_path), "dgCMatrix")
genes <- read.delim(genes_path, stringsAsFactors = FALSE, check.names = FALSE)[[1]]
if (length(genes) != nrow(counts)) stop("Gene index length does not match matrix rows")
rownames(counts) <- make.unique(as.character(genes))
if (!target_tf %in% rownames(counts)) stop(target_tf, " is absent from the network gene matrix")
nc_ncells_used <- min(nc_ncells_requested, ncol(counts))
if (nc_ncells_used < 1L) stop("No cells are available")

extract_diff_regulation <- function(result, tf) {
  candidates <- c("diffRegulation", "diff_regulation", "diffRegulationResult")
  table <- NULL
  for (name in candidates) {
    if (!is.null(result[[name]]) && is.data.frame(result[[name]])) {
      table <- as.data.frame(result[[name]])
      break
    }
  }
  if (is.null(table)) stop("scTenifoldKnk result lacks differential-regulation table")
  if (!"gene" %in% names(table)) table$gene <- rownames(table)
  if (!"tf" %in% names(table)) table$tf <- tf
  aliases <- list(
    distance = c("distance", "Distance"),
    p.adj = c("p.adj", "p_adj", "FDR"),
    p.value = c("p.value", "pvalue", "P.Value"),
    Z = c("Z", "z", "zscore"),
    FC = c("FC", "fc", "fold_change")
  )
  for (canonical in names(aliases)) {
    if (canonical %in% names(table)) next
    hit <- aliases[[canonical]][aliases[[canonical]] %in% names(table)]
    if (length(hit)) table[[canonical]] <- table[[hit[1]]]
  }
  required <- c("tf", "gene", "distance", "p.adj", "p.value", "Z", "FC")
  missing <- setdiff(required, names(table))
  if (length(missing)) stop("Missing result columns: ", paste(missing, collapse = ", "))
  table[, required, drop = FALSE]
}

network_to_edge_summary <- function(network, tf, network_label) {
  if (is.null(network)) return(data.frame())
  if (inherits(network, "Matrix") || inherits(network, "dgCMatrix") || inherits(network, "dgTMatrix") || inherits(network, "matrix")) {
    if (inherits(network, "matrix")) {
      values <- as.numeric(network)
      n_edges <- sum(is.finite(values) & values != 0)
      mean_abs_weight <- if (n_edges) mean(abs(values[is.finite(values) & values != 0])) else 0
      max_abs_weight <- if (n_edges) max(abs(values[is.finite(values) & values != 0])) else 0
    } else {
      edges <- Matrix::summary(as(network, "dgCMatrix"))
      n_edges <- nrow(edges)
      mean_abs_weight <- if (n_edges) mean(abs(edges$x)) else 0
      max_abs_weight <- if (n_edges) max(abs(edges$x)) else 0
    }
    row_names <- rownames(network)
    col_names <- colnames(network)
    n_nodes <- length(unique(c(row_names, col_names)))
    if (!n_nodes) n_nodes <- max(nrow(network), ncol(network))
    return(data.frame(
      tf = tf,
      network = network_label,
      n_nodes = n_nodes,
      n_edges = n_edges,
      mean_abs_weight = mean_abs_weight,
      max_abs_weight = max_abs_weight,
      stringsAsFactors = FALSE
    ))
  }
  if (is.data.frame(network)) {
    weight_col <- intersect(c("weight", "Weight", "value", "coef", "importance"), names(network))
    weights <- if (length(weight_col)) as.numeric(network[[weight_col[1]]]) else rep(1, nrow(network))
    node_cols <- intersect(c("source", "target", "from", "to", "gene1", "gene2"), names(network))
    nodes <- if (length(node_cols)) unique(unlist(network[node_cols], use.names = FALSE)) else character()
    return(data.frame(
      tf = tf,
      network = network_label,
      n_nodes = length(nodes),
      n_edges = nrow(network),
      mean_abs_weight = if (length(weights)) mean(abs(weights), na.rm = TRUE) else 0,
      max_abs_weight = if (length(weights)) max(abs(weights), na.rm = TRUE) else 0,
      stringsAsFactors = FALSE
    ))
  }
  data.frame()
}

extract_network_summaries <- function(result, tf) {
  if (is.null(result$tensorNetworks)) {
    return(data.frame(tf = tf, network = NA_character_, n_nodes = NA_integer_, n_edges = NA_integer_,
                      mean_abs_weight = NA_real_, max_abs_weight = NA_real_, stringsAsFactors = FALSE))
  }
  rows <- list()
  for (name in names(result$tensorNetworks)) {
    summary <- network_to_edge_summary(result$tensorNetworks[[name]], tf, paste0("tensorNetworks.", name))
    if (nrow(summary)) rows[[name]] <- summary
  }
  if (!length(rows)) {
    return(data.frame(tf = tf, network = NA_character_, n_nodes = NA_integer_, n_edges = NA_integer_,
                      mean_abs_weight = NA_real_, max_abs_weight = NA_real_, stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

started <- Sys.time()
set.seed(seed)
result <- scTenifoldKnk::scTenifoldKnk(
  countMatrix = counts,
  gKO = target_tf,
  qc = qc,
  qc_minCells = qc_min_cells,
  nc_nNet = nc_nnet,
  nc_nCells = nc_ncells_used,
  nc_nComp = nc_ncomp,
  ma_nDim = ma_ndim,
  nCores = ncores
)

perturbation <- extract_diff_regulation(result, target_tf)
perturbation$tf <- target_tf
perturbation$subset <- subset
perturbation$axis <- axis
perturbation$seed <- seed
network_summary <- extract_network_summaries(result, target_tf)
network_summary$subset <- subset
network_summary$axis <- axis
network_summary$seed <- seed
result_names <- names(result)
contract <- data.frame(
  tf = target_tf,
  axis = axis,
  subset = subset,
  seed = seed,
  has_wt_network = !is.null(result$tensorNetworks$WT),
  has_ko_network = !is.null(result$tensorNetworks$KO),
  has_manifold_alignment = !is.null(result$manifoldAlignment),
  has_perturbation_table = nrow(perturbation) > 0,
  n_perturbation_genes = nrow(perturbation),
  result_names = paste(result_names, collapse = ";"),
  stringsAsFactors = FALSE
)

prefix <- paste0(target_tf, "_", subset, "_seed", seed)
rds_path <- file.path(output_dir, paste0(prefix, "_result.rds"))
perturbation_path <- file.path(metadata_dir, paste0(prefix, "_perturbation_genes.tsv"))
network_path <- file.path(metadata_dir, paste0(prefix, "_network_adjacency_summary.tsv"))
contract_path <- file.path(metadata_dir, paste0(prefix, "_result_contract.tsv"))
report_path <- file.path(metadata_dir, paste0(prefix, "_run_report.json"))
saveRDS(result, rds_path)
write.table(perturbation, perturbation_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(network_summary, network_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contract, contract_path, sep = "\t", quote = FALSE, row.names = FALSE)

report <- list(
  module = "scTenifoldKnk reproducibility audit v2",
  target_tf = target_tf,
  axis = axis,
  subset = subset,
  seed = seed,
  status = "complete",
  matrix_orientation = "genes_x_cells",
  input = list(matrix = path_for_report(matrix_path), genes = path_for_report(genes_path)),
  n_genes = nrow(counts),
  n_cells_total = ncol(counts),
  n_cells_used = nc_ncells_used,
  parameters = list(qc = qc, qc_minCells = qc_min_cells, nc_nNet = nc_nnet,
                    nc_nCells_requested = nc_ncells_requested, nc_nCells_used = nc_ncells_used,
                    nc_nComp = nc_ncomp, ma_nDim = ma_ndim, nCores = ncores, seed = seed),
  n_significant_excluding_target = sum(perturbation$gene != target_tf &
    is.finite(as.numeric(perturbation$p.adj)) & as.numeric(perturbation$p.adj) < 0.05),
  result_structure = list(names = result_names,
                          has_wt_network = contract$has_wt_network,
                          has_ko_network = contract$has_ko_network,
                          has_manifold_alignment = contract$has_manifold_alignment,
                          has_perturbation_table = contract$has_perturbation_table),
  runtime = list(r = R.version.string,
                 scTenifoldKnk = as.character(packageVersion("scTenifoldKnk")),
                 Matrix = as.character(packageVersion("Matrix")),
                 jsonlite = as.character(packageVersion("jsonlite"))),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  outputs = list(result_rds = path_for_report(rds_path),
                 perturbation_genes = path_for_report(perturbation_path),
                 network_adjacency_summary = path_for_report(network_path),
                 result_contract = path_for_report(contract_path),
                 report = path_for_report(report_path)),
  caveat = "scTenifoldKnk is computational virtual network perturbation evidence; manifold distance is unsigned and is not interpreted as biological activation or suppression."
)
write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
message("Completed ", target_tf, " seed ", seed, "; report=", report_path)
