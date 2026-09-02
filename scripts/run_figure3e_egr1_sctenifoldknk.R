#!/usr/bin/env Rscript

## Run EGR1 scTenifoldKnk with a dedicated Figure 3 output contract.

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    input_dir = "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative",
    output_dir = "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/results",
    metadata_dir = "metadata/driver/figure3e_egr1",
    subset = "stressed_regenerative",
    seeds = "15071990,15071991,15071992",
    nc_nnet = 10,
    nc_ncells = 500,
    nc_ncomp = 3,
    ma_ndim = 2,
    ncores = 8,
    qc = "false",
    qc_min_cells = 3
  )
  for (item in args) {
    parts <- strsplit(sub("^--", "", item), "=", fixed = TRUE)[[1]]
    if (length(parts) == 2L && parts[1] %in% names(out)) out[[parts[1]]] <- parts[2]
  }
  out$seeds <- as.integer(trimws(strsplit(out$seeds, ",", fixed = TRUE)[[1]]))
  out$nc_nnet <- as.integer(out$nc_nnet)
  out$nc_ncells <- as.integer(out$nc_ncells)
  out$nc_ncomp <- as.integer(out$nc_ncomp)
  out$ma_ndim <- as.integer(out$ma_ndim)
  out$ncores <- as.integer(out$ncores)
  out$qc <- tolower(as.character(out$qc)) %in% c("true", "t", "1", "yes", "y")
  out$qc_min_cells <- as.integer(out$qc_min_cells)
  out
}

extract_diff_regulation <- function(result, tf = "EGR1") {
  candidates <- c("diffRegulation", "diff_regulation", "diffRegulationResult")
  table <- NULL
  for (name in candidates) {
    if (!is.null(result[[name]]) && is.data.frame(result[[name]])) {
      table <- as.data.frame(result[[name]])
      break
    }
  }
  if (is.null(table)) stop("scTenifoldKnk result lacks a differential-regulation table")
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
  if (length(missing)) stop("Missing scTenifoldKnk result columns: ", paste(missing, collapse = ", "))
  table[, required, drop = FALSE]
}

consensus_table <- function(combined, expected_seeds) {
  split_rows <- split(combined, combined$gene)
  rows <- lapply(split_rows, function(frame) {
    frame$distance <- as.numeric(frame$distance)
    frame$p.adj <- as.numeric(frame$p.adj)
    frame$p.value <- as.numeric(frame$p.value)
    frame$Z <- as.numeric(frame$Z)
    frame$FC <- as.numeric(frame$FC)
    data.frame(
      tf = "EGR1",
      gene = as.character(frame$gene[1]),
      distance = median(frame$distance, na.rm = TRUE),
      distance_min = min(frame$distance, na.rm = TRUE),
      distance_max = max(frame$distance, na.rm = TRUE),
      p.adj = max(frame$p.adj, na.rm = TRUE),
      p.value = max(frame$p.value, na.rm = TRUE),
      Z = median(frame$Z, na.rm = TRUE),
      FC = median(frame$FC, na.rm = TRUE),
      n_successful_seeds = length(unique(frame$seed)),
      n_expected_seeds = expected_seeds,
      significant_seed_fraction = mean(frame$p.adj < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$distance, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

main <- function() {
  args <- parse_args()
  if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
    stop("R package scTenifoldKnk is required for Figure 3E")
  }
  input_dir <- normalizePath(args$input_dir, mustWork = TRUE)
  output_dir <- args$output_dir
  metadata_dir <- args$metadata_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_path <- file.path(
    input_dir,
    paste0("figure3e_egr1_", args$subset, "_counts_genes_x_cells.mtx")
  )
  genes_path <- file.path(input_dir, paste0("figure3e_egr1_", args$subset, "_genes.tsv"))
  if (!all(file.exists(c(matrix_path, genes_path)))) stop("Figure 3E matrix or gene index is missing")
  counts <- Matrix::readMM(matrix_path)
  genes <- read.delim(genes_path, stringsAsFactors = FALSE)[[1]]
  rownames(counts) <- make.unique(as.character(genes))
  counts <- as(counts, "dgCMatrix")
  if (!"EGR1" %in% rownames(counts)) stop("EGR1 is absent from the network gene matrix")
  actual_ncells <- min(args$nc_ncells, ncol(counts))

  seed_tables <- list()
  run_rows <- list()
  failures <- list()
  for (seed in args$seeds) {
    started <- Sys.time()
    message("Running Figure 3E EGR1 scTenifoldKnk: subset=", args$subset, ", seed=", seed)
    set.seed(seed)
    result <- tryCatch(
      scTenifoldKnk::scTenifoldKnk(
        countMatrix = counts,
        gKO = "EGR1",
        qc = args$qc,
        qc_minCells = args$qc_min_cells,
        nc_nNet = args$nc_nnet,
        nc_nCells = actual_ncells,
        nc_nComp = args$nc_ncomp,
        ma_nDim = args$ma_ndim,
        nCores = args$ncores
      ),
      error = function(error) error
    )
    finished <- Sys.time()
    if (inherits(result, "error")) {
      failures[[as.character(seed)]] <- conditionMessage(result)
      run_rows[[as.character(seed)]] <- data.frame(
        subset = args$subset,
        seed = seed,
        status = "failed",
        elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
        result_rds = "",
        perturbation_tsv = "",
        error = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      next
    }
    table <- extract_diff_regulation(result, "EGR1")
    table$subset <- args$subset
    table$seed <- seed
    stem <- paste0("figure3e_egr1_", args$subset, "_seed", seed)
    rds_path <- file.path(output_dir, paste0(stem, "_result.rds"))
    tsv_path <- file.path(output_dir, paste0(stem, "_perturbation_genes.tsv"))
    saveRDS(result, rds_path)
    write.table(table, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
    seed_tables[[as.character(seed)]] <- table
    run_rows[[as.character(seed)]] <- data.frame(
      subset = args$subset,
      seed = seed,
      status = "success",
      elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
      result_rds = normalizePath(rds_path, mustWork = FALSE),
      perturbation_tsv = normalizePath(tsv_path, mustWork = FALSE),
      error = "",
      stringsAsFactors = FALSE
    )
  }
  run_log <- do.call(rbind, run_rows)
  run_log_path <- file.path(metadata_dir, paste0("figure3e_egr1_", args$subset, "_run_log.tsv"))
  write.table(run_log, run_log_path, sep = "\t", quote = FALSE, row.names = FALSE)

  if (!length(seed_tables)) {
    report <- list(
      module = "Figure 3E scTenifoldKnk",
      target_tf = "EGR1",
      subset = args$subset,
      status = "failed",
      parameters = list(
        nc_nNet = args$nc_nnet,
        nc_nCells_requested = args$nc_ncells,
        nc_nCells_used = actual_ncells,
        seeds = args$seeds,
        nCores = args$ncores
      ),
      failures = failures,
      run_log = normalizePath(run_log_path, mustWork = FALSE)
    )
    report_path <- file.path(metadata_dir, paste0("figure3e_egr1_", args$subset, "_run_report.json"))
    jsonlite::write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE)
    stop("All EGR1 scTenifoldKnk seeds failed for subset ", args$subset)
  }

  combined <- do.call(rbind, seed_tables)
  combined_path <- file.path(
    metadata_dir,
    paste0("figure3e_egr1_", args$subset, "_all_seed_perturbation_genes.tsv")
  )
  write.table(combined, combined_path, sep = "\t", quote = FALSE, row.names = FALSE)
  consensus <- consensus_table(combined, length(args$seeds))
  consensus$subset <- args$subset
  consensus_path <- file.path(
    metadata_dir,
    paste0("figure3e_egr1_", args$subset, "_consensus_perturbation_genes.tsv")
  )
  write.table(consensus, consensus_path, sep = "\t", quote = FALSE, row.names = FALSE)

  review_risks <- list()
  if (args$nc_nnet < 10L) {
    review_risks[[length(review_risks) + 1L]] <- list(
      flag = "nc_nNet_below_formal_contract",
      severity = "review_attention",
      detail = paste0("nc_nNet=", args$nc_nnet, " was used.")
    )
  }
  if (actual_ncells < 500L) {
    review_risks[[length(review_risks) + 1L]] <- list(
      flag = "nc_nCells_below_formal_contract",
      severity = "review_attention",
      detail = paste0("nc_nCells=", actual_ncells, " was limited by subset size.")
    )
  }
  if (length(seed_tables) < 2L) {
    review_risks[[length(review_risks) + 1L]] <- list(
      flag = "multiple_seed_replication_not_met",
      severity = "review_attention",
      detail = paste0("Only ", length(seed_tables), " successful seed(s).")
    )
  }
  report <- list(
    module = "Figure 3E scTenifoldKnk",
    method = "EGR1 virtual knockout with conservative across-seed consensus",
    target_tf = "EGR1",
    subset = args$subset,
    status = if (length(seed_tables) == length(args$seeds)) "complete" else "partial",
    matrix_orientation = "genes_x_cells",
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    parameters = list(
      qc = args$qc,
      qc_minCells = args$qc_min_cells,
      nc_nNet = args$nc_nnet,
      nc_nCells_requested = args$nc_ncells,
      nc_nCells_used = actual_ncells,
      nc_nComp = args$nc_ncomp,
      ma_nDim = args$ma_ndim,
      nCores = args$ncores,
      seeds = args$seeds
    ),
    n_successful_seeds = length(seed_tables),
    n_failed_seeds = length(failures),
    consensus = list(
      distance = "median across successful seeds",
      p_adjust = "maximum p.adj across successful seeds",
      p_value = "maximum p.value across successful seeds",
      significance = "conservative p.adj < 0.05 after maximum-across-seed aggregation",
      n_significant_excluding_egr1 = sum(consensus$gene != "EGR1" & consensus$p.adj < 0.05, na.rm = TRUE)
    ),
    review_risk_flags = review_risks,
    failures = failures,
    outputs = list(
      all_seed_results = normalizePath(combined_path, mustWork = FALSE),
      consensus_results = normalizePath(consensus_path, mustWork = FALSE),
      run_log = normalizePath(run_log_path, mustWork = FALSE)
    ),
    caveat = "scTenifoldKnk is a computational virtual knockout. Distance is an unsigned network-displacement statistic and is not interpreted as up- or down-regulation."
  )
  report_path <- file.path(metadata_dir, paste0("figure3e_egr1_", args$subset, "_run_report.json"))
  report$outputs$report <- normalizePath(report_path, mustWork = FALSE)
  jsonlite::write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  message("Figure 3E scTenifoldKnk report written: ", report_path)
}

main()

