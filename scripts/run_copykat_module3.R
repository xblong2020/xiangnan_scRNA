suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

parse_args <- function(argv) {
  args <- list(
    manifest = file.path("metadata", "copykat_module3", "copykat_input_manifest.tsv"),
    metadata_dir = file.path("metadata", "copykat_module3"),
    runs = "",
    overwrite = FALSE,
    genome = "hg20",
    n_cores = 1L,
    ks_cut = 0.1,
    dense_retry_max_cells = 6000L,
    sparse_first = FALSE,
    skip_heatmap = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key == "--overwrite") {
      args$overwrite <- TRUE
      i <- i + 1L
      next
    }
    if (key == "--sparse-first") {
      args$sparse_first <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop("Missing value for ", key)
    val <- argv[[i + 1L]]
    if (key == "--manifest") args$manifest <- val
    else if (key == "--metadata-dir") args$metadata_dir <- val
    else if (key == "--runs") args$runs <- val
    else if (key == "--genome") args$genome <- val
    else if (key == "--n-cores") args$n_cores <- as.integer(val)
    else if (key == "--ks-cut") args$ks_cut <- as.numeric(val)
    else if (key == "--dense-retry-max-cells") args$dense_retry_max_cells <- as.integer(val)
    else stop("Unknown argument: ", key)
    i <- i + 2L
  }
  args
}

resolve_path <- function(path) {
  path <- as.character(path)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  candidate <- file.path(getwd(), path)
  if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_lines <- function(path) {
  if (!file.exists(path)) return(character())
  trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
}

read_run_matrix <- function(run_dir) {
  matrix_path <- file.path(run_dir, "raw_counts_gene_by_cell.mtx")
  genes_path <- file.path(run_dir, "genes.tsv")
  cell_map_path <- file.path(run_dir, "cell_map.tsv")
  if (!file.exists(matrix_path)) stop("Missing matrix: ", matrix_path)
  if (!file.exists(genes_path)) stop("Missing genes: ", genes_path)
  if (!file.exists(cell_map_path)) stop("Missing cell map: ", cell_map_path)

  rawmat <- Matrix::readMM(matrix_path)
  rawmat <- as(rawmat, "CsparseMatrix")
  genes <- read.delim(genes_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")$gene
  cell_map <- read.delim(cell_map_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  if (length(genes) != nrow(rawmat)) stop("Gene count does not match matrix rows in ", run_dir)
  if (nrow(cell_map) != ncol(rawmat)) stop("Cell count does not match matrix columns in ", run_dir)
  rownames(rawmat) <- make.unique(as.character(genes))
  colnames(rawmat) <- as.character(cell_map$cell_key)
  list(rawmat = rawmat, cell_map = cell_map)
}

make_copykat_no_heatmap <- function() {
  fn <- copykat::copykat
  skip_expr <- quote({
    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
    names(reslts) <- c("CNAmat", "hclustering")
    return(reslts)
  })
  rewrite <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("{"))) {
      new_expr <- list(as.name("{"))
      if (length(expr) >= 2L) {
        for (j in 2L:length(expr)) {
          child <- expr[[j]]
          if (
            is.call(child) &&
              identical(child[[1L]], as.name("print")) &&
              length(child) >= 2L &&
              is.character(child[[2L]]) &&
              grepl("ploting heatmap", child[[2L]], fixed = TRUE)
          ) {
            return(as.call(c(new_expr, as.list(skip_expr)[-1L])))
          }
          new_expr[[length(new_expr) + 1L]] <- rewrite(child)
        }
      }
      return(as.call(new_expr))
    }
    if (is.call(expr) && length(expr) >= 2L) {
      for (j in 2L:length(expr)) {
        expr[[j]] <- rewrite(expr[[j]])
      }
    }
    expr
  }
  body(fn) <- rewrite(body(fn))
  fn
}

run_copykat_once <- function(rawmat, norm_cell_names, run_id, args, use_dense) {
  mat <- rawmat
  if (use_dense) {
    mat <- as.matrix(rawmat)
    storage.mode(mat) <- "numeric"
  }
  copykat_fn <- if (args$skip_heatmap) make_copykat_no_heatmap() else copykat::copykat
  copykat_fn(
    rawmat = mat,
    id.type = "S",
    cell.line = "no",
    ngene.chr = 5,
    min.gene.per.cell = 200,
    LOW.DR = 0.05,
    UP.DR = 0.1,
    win.size = 25,
    norm.cell.names = norm_cell_names,
    KS.cut = args$ks_cut,
    sam.name = run_id,
    distance = "euclidean",
    output.seg = "FALSE",
    plot.genes = "FALSE",
    genome = args$genome,
    n.cores = args$n_cores
  )
}

run_one <- function(row, args) {
  run_id <- as.character(row$run_id)
  run_dir <- resolve_path(row$run_dir)
  prediction_path <- file.path(run_dir, paste0(run_id, "_copykat_prediction.txt"))
  result_path <- file.path(run_dir, paste0(run_id, "_copykat_result.rds"))
  log_path <- file.path(run_dir, paste0(run_id, "_copykat_runner.log"))
  if (file.exists(prediction_path) && !args$overwrite) {
    return(data.frame(
      run_id = run_id, status = if (file.exists(result_path)) "exists" else "exists_prediction_only", run_dir = run_dir,
      prediction_path = prediction_path, result_rds = if (file.exists(result_path)) result_path else "",
      n_cells = as.integer(row$n_cells %||% NA), n_candidate = as.integer(row$n_candidate %||% NA),
      n_reference = as.integer(row$n_reference %||% NA), message = "",
      stringsAsFactors = FALSE
    ))
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  message("COPYKAT ", run_id, " cells=", row$n_cells, " candidates=", row$n_candidate, " refs=", row$n_reference)
  loaded <- read_run_matrix(run_dir)
  norm_cell_names <- read_lines(file.path(run_dir, "normal_cell_keys.txt"))
  norm_cell_names <- norm_cell_names[norm_cell_names %in% colnames(loaded$rawmat)]
  if (length(norm_cell_names) < 2L) stop(run_id, ": CopyKAT requires at least two normal reference cells")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run_dir)
  out_con <- file(log_path, open = "w", encoding = "UTF-8")
  msg_con <- file(log_path, open = "a", encoding = "UTF-8")
  sink(out_con, split = TRUE)
  sink(msg_con, type = "message")
  on.exit({
    if (sink.number(type = "message") > 0L) sink(type = "message")
    if (sink.number() > 0L) sink()
    close(msg_con)
    close(out_con)
  }, add = TRUE)

  used_dense <- !args$sparse_first
  if (args$sparse_first) {
    result <- tryCatch(
      run_copykat_once(loaded$rawmat, norm_cell_names, run_id, args, use_dense = FALSE),
      error = function(err) {
        if (ncol(loaded$rawmat) <= args$dense_retry_max_cells) {
          message("Sparse CopyKAT failed for ", run_id, "; retrying dense. Error: ", conditionMessage(err))
          used_dense <<- TRUE
          run_copykat_once(loaded$rawmat, norm_cell_names, run_id, args, use_dense = TRUE)
        } else {
          stop(err)
        }
      }
    )
  } else {
    result <- run_copykat_once(loaded$rawmat, norm_cell_names, run_id, args, use_dense = TRUE)
  }
  saveRDS(result, result_path)
  if (!file.exists(prediction_path)) {
    candidates <- list.files(run_dir, pattern = "copykat_prediction\\.txt$", full.names = TRUE)
    if (length(candidates) > 0L) prediction_path <- candidates[[1L]]
  }
  data.frame(
    run_id = run_id,
    status = if (file.exists(prediction_path)) "ok" else "missing_prediction",
    run_dir = run_dir,
    prediction_path = if (file.exists(prediction_path)) prediction_path else "",
    result_rds = result_path,
    n_cells = ncol(loaded$rawmat),
    n_candidate = sum(loaded$cell_map$cnv_role == "candidate"),
    n_reference = sum(loaded$cell_map$cnv_role == "reference"),
    message = paste0("used_dense=", used_dense),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  args$manifest <- resolve_path(args$manifest)
  args$metadata_dir <- resolve_path(args$metadata_dir)
  dir.create(args$metadata_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- read.delim(args$manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  if (nrow(manifest) == 0L) stop("Manifest has no runs: ", args$manifest)
  if (nzchar(args$runs)) {
    keep <- trimws(strsplit(args$runs, ",", fixed = TRUE)[[1L]])
    manifest <- manifest[manifest$run_id %in% keep, , drop = FALSE]
  }
  if (nrow(manifest) == 0L) stop("No requested runs found in manifest")

  status_rows <- list()
  status_path <- file.path(args$metadata_dir, "copykat_run_status.tsv")
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    status_rows[[length(status_rows) + 1L]] <- tryCatch(
      run_one(row, args),
      error = function(err) {
        data.frame(
          run_id = as.character(row$run_id),
          status = "error",
          run_dir = as.character(row$run_dir),
          prediction_path = "",
          result_rds = "",
          n_cells = suppressWarnings(as.integer(row$n_cells %||% NA)),
          n_candidate = suppressWarnings(as.integer(row$n_candidate %||% NA)),
          n_reference = suppressWarnings(as.integer(row$n_reference %||% NA)),
          message = conditionMessage(err),
          stringsAsFactors = FALSE
        )
      }
    )
    status <- do.call(rbind, status_rows)
    write.table(status, status_path, sep = "\t", row.names = FALSE, quote = FALSE)
    gc(verbose = FALSE)
  }
  status <- do.call(rbind, status_rows)
  write.table(status, status_path, sep = "\t", row.names = FALSE, quote = FALSE)
  print(status)
  if (any(status$status %in% c("error", "missing_prediction"))) {
    quit(status = 1L)
  }
  invisible(status)
}

main()
