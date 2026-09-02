suppressPackageStartupMessages({
  library(Matrix)
})

default_tfs <- c(
  "ATF3", "CEBPB", "EGR1", "FOS", "HLF", "HNF4A", "IRF1", "JUN",
  "JUNB", "JUND", "MAFB", "MAFF", "MYC", "PPARA", "SOX4"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    input_dir = "data/processed/driver/sctenifoldknk_module7_1/driver_union_all",
    output_dir = "data/processed/driver/sctenifoldknk_module7_2/driver_union_all",
    metadata_dir = "metadata/driver",
    subset = "driver_union_all",
    matrix_file = "sctenifoldknk_counts_genes_x_cells.mtx",
    genes_file = "sctenifoldknk_genes.tsv",
    tfs = paste(default_tfs, collapse = ","),
    qc = "false",
    qc_min_cells = 3,
    nc_nnet = 10,
    nc_ncells = 500,
    nc_ncomp = 3,
    ma_ndim = 2,
    ncores = 1,
    seed = 1
  )
  if (length(args) == 0) {
    return(out)
  }
  for (item in args) {
    parts <- strsplit(sub("^--", "", item), "=", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] %in% names(out)) {
      out[[parts[1]]] <- parts[2]
    }
  }
  out$qc_min_cells <- as.integer(out$qc_min_cells)
  out$qc <- tolower(as.character(out$qc)) %in% c("true", "t", "1", "yes", "y")
  out$nc_nnet <- as.integer(out$nc_nnet)
  out$nc_ncells <- as.integer(out$nc_ncells)
  out$nc_ncomp <- as.integer(out$nc_ncomp)
  out$ma_ndim <- as.integer(out$ma_ndim)
  out$ncores <- as.integer(out$ncores)
  out$seed <- as.integer(out$seed)
  out$tfs <- trimws(strsplit(out$tfs, ",")[[1]])
  out
}

extract_diff_regulation <- function(result, tf) {
  candidates <- c("diffRegulation", "diff_regulation", "diffRegulationResult")
  table <- NULL
  for (name in candidates) {
    if (!is.null(result[[name]]) && is.data.frame(result[[name]])) {
      table <- result[[name]]
      break
    }
  }
  if (is.null(table)) {
    stop("scTenifoldKnk result does not contain a diffRegulation data.frame for ", tf)
  }
  table <- as.data.frame(table)
  if (!"gene" %in% colnames(table)) {
    table$gene <- rownames(table)
  }
  if (!"tf" %in% colnames(table)) {
    table$tf <- tf
  }
  preferred <- c("tf", "gene", "distance", "Distance", "p.adj", "p_adj", "pvalue", "p.value", "FDR")
  table[, unique(c(intersect(preferred, colnames(table)), setdiff(colnames(table), preferred))), drop = FALSE]
}

network_to_edge_summary <- function(network, tf, network_label) {
  if (is.null(network)) {
    return(data.frame())
  }
  if (inherits(network, "dgCMatrix") || inherits(network, "matrix")) {
    mat <- as(network, "dgCMatrix")
    edges <- Matrix::summary(mat)
    row_names <- rownames(mat)
    col_names <- colnames(mat)
    if (is.null(row_names)) {
      row_names <- as.character(seq_len(nrow(mat)))
    }
    if (is.null(col_names)) {
      col_names <- as.character(seq_len(ncol(mat)))
    }
    return(data.frame(
      tf = tf,
      network = network_label,
      n_nodes = length(unique(c(row_names, col_names))),
      n_edges = nrow(edges),
      mean_abs_weight = ifelse(nrow(edges) > 0, mean(abs(edges$x)), 0),
      max_abs_weight = ifelse(nrow(edges) > 0, max(abs(edges$x)), 0),
      stringsAsFactors = FALSE
    ))
  }
  if (is.data.frame(network)) {
    weight_col <- intersect(c("weight", "Weight", "value", "coef", "importance"), colnames(network))
    weights <- if (length(weight_col) > 0) as.numeric(network[[weight_col[[1]]]]) else rep(1, nrow(network))
    node_cols <- intersect(c("source", "target", "from", "to", "gene1", "gene2"), colnames(network))
    nodes <- unique(unlist(network[node_cols], use.names = FALSE))
    return(data.frame(
      tf = tf,
      network = network_label,
      n_nodes = length(nodes),
      n_edges = nrow(network),
      mean_abs_weight = ifelse(length(weights) > 0, mean(abs(weights), na.rm = TRUE), 0),
      max_abs_weight = ifelse(length(weights) > 0, max(abs(weights), na.rm = TRUE), 0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame()
}

extract_network_summaries <- function(result, tf) {
  network_names <- names(result)[grepl("net|network|wt|wild|ko|knock", names(result), ignore.case = TRUE)]
  rows <- list()
  for (name in network_names) {
    item <- result[[name]]
    if (is.list(item) && !is.data.frame(item)) {
      for (child_name in names(item)) {
        child_summary <- network_to_edge_summary(item[[child_name]], tf, paste(name, child_name, sep = "."))
        if (nrow(child_summary) > 0) {
          rows[[paste(name, child_name, sep = ".")]] <- child_summary
        }
      }
    } else {
      summary <- network_to_edge_summary(item, tf, name)
      if (nrow(summary) > 0) {
        rows[[name]] <- summary
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(
      tf = tf,
      network = NA_character_,
      n_nodes = NA_integer_,
      n_edges = NA_integer_,
      mean_abs_weight = NA_real_,
      max_abs_weight = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

summarize_result_contract <- function(result, tf, perturbation_table) {
  names_present <- names(result)
  data.frame(
    tf = tf,
    has_wt_network = any(grepl("wt|wild|net1|network", names_present, ignore.case = TRUE)),
    has_ko_network = any(grepl("ko|knock|net2|network", names_present, ignore.case = TRUE)),
    has_manifold_alignment = any(grepl("manifold|alignment|ma_", names_present, ignore.case = TRUE)),
    has_perturbation_table = is.data.frame(perturbation_table) && nrow(perturbation_table) > 0,
    n_perturbation_genes = nrow(perturbation_table),
    result_names = paste(names_present, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args()
  if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
    stop("R package scTenifoldKnk is required for Module 7.2")
  }

  set.seed(args$seed)
  input_dir <- normalizePath(args$input_dir, mustWork = TRUE)
  output_dir <- args$output_dir
  metadata_dir <- args$metadata_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_path <- file.path(input_dir, args$matrix_file)
  genes_path <- file.path(input_dir, args$genes_file)
  count_matrix <- Matrix::readMM(matrix_path)
  genes <- read.delim(genes_path, stringsAsFactors = FALSE)[[1]]
  rownames(count_matrix) <- make.unique(as.character(genes))
  count_matrix <- as(count_matrix, "dgCMatrix")

  missing_tfs <- setdiff(args$tfs, rownames(count_matrix))
  if (length(missing_tfs) > 0) {
    stop("Missing knockout TFs in count matrix: ", paste(missing_tfs, collapse = ", "))
  }

  all_tables <- list()
  contract_rows <- list()
  network_summary_rows <- list()
  log_rows <- list()
  for (tf in args$tfs) {
    started <- Sys.time()
    message("Running scTenifoldKnk knockout: ", tf)
    result <- scTenifoldKnk::scTenifoldKnk(
      countMatrix = count_matrix,
      gKO = tf,
      qc = args$qc,
      qc_minCells = args$qc_min_cells,
      nc_nNet = args$nc_nnet,
      nc_nCells = args$nc_ncells,
      nc_nComp = args$nc_ncomp,
      ma_nDim = args$ma_ndim,
      nCores = args$ncores
    )
    rds_path <- file.path(output_dir, paste0("sctenifoldknk_", args$subset, "_", tf, ".rds"))
    saveRDS(result, rds_path)
    perturbation <- extract_diff_regulation(result, tf)
    perturbation$subset <- args$subset
    tsv_path <- file.path(output_dir, paste0("sctenifoldknk_", args$subset, "_", tf, "_perturbation_genes.tsv"))
    write.table(perturbation, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
    network_summary <- extract_network_summaries(result, tf)
    network_summary$subset <- args$subset
    network_summary_path <- file.path(output_dir, paste0("sctenifoldknk_", args$subset, "_", tf, "_network_adjacency_summary.tsv"))
    write.table(network_summary, network_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
    all_tables[[tf]] <- perturbation
    contract_rows[[tf]] <- summarize_result_contract(result, tf, perturbation)
    network_summary_rows[[tf]] <- network_summary
    log_rows[[tf]] <- data.frame(
      tf = tf,
      subset = args$subset,
      started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
      finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      result_rds = rds_path,
      perturbation_tsv = tsv_path,
      network_adjacency_summary_tsv = network_summary_path,
      stringsAsFactors = FALSE
    )
  }

  combined <- do.call(rbind, all_tables)
  contract <- do.call(rbind, contract_rows)
  network_summary_all <- do.call(rbind, network_summary_rows)
  run_log <- do.call(rbind, log_rows)
  combined_path <- file.path(metadata_dir, paste0("sctenifoldknk_module7_2_", args$subset, "_perturbation_genes.tsv"))
  contract_path <- file.path(metadata_dir, paste0("sctenifoldknk_module7_2_", args$subset, "_result_contract.tsv"))
  network_summary_path <- file.path(metadata_dir, paste0("sctenifoldknk_module7_2_", args$subset, "_network_adjacency_summary.tsv"))
  log_path <- file.path(metadata_dir, paste0("sctenifoldknk_module7_2_", args$subset, "_run_log.tsv"))
  report_path <- file.path(metadata_dir, paste0("sctenifoldknk_module7_2_", args$subset, "_report.json"))
  write.table(combined, combined_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(contract, contract_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(network_summary_all, network_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(run_log, log_path, sep = "\t", quote = FALSE, row.names = FALSE)

  report <- list(
    module = "7.2",
    method = "Batch scTenifoldKnk virtual knockout for CellOracle TFs",
    subset = args$subset,
    matrix_orientation = "genes_x_cells",
    n_genes = nrow(count_matrix),
    n_cells = ncol(count_matrix),
    tfs = args$tfs,
    parameters = list(
      qc_minCells = args$qc_min_cells,
      qc = args$qc,
      nc_nNet = args$nc_nnet,
      nc_nCells = args$nc_ncells,
      nc_nComp = args$nc_ncomp,
      ma_nDim = args$ma_ndim,
      nCores = args$ncores,
      seed = args$seed
    ),
    outputs = list(
      perturbation_genes = combined_path,
      result_contract = contract_path,
      network_adjacency_summary = network_summary_path,
      run_log = log_path,
      report = report_path
    )
  )
  jsonlite::write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE)
  message("Wrote Module 7.2 report: ", report_path)
}

main()
