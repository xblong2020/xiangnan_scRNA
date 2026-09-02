args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(doParallel)
  library(dplyr)
  library(HiClimR)
  library(magrittr)
  library(Matrix)
  library(plyr)
  library(Rfast)
  library(RSpectra)
  library(stringr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/run_figure1c_cytotrace2.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    input_manifest = file.path(ROOT, "metadata/figure1c/figure1c_cytotrace2_input_manifest.tsv"),
    output_dir = file.path(ROOT, "metadata/figure1c"),
    source_cache = file.path(ROOT, "tmp/cytotrace2_r_source"),
    mode = "all",
    species = "human",
    batch_size = 5000,
    smooth_batch_size = 1000,
    seed = 20260708
  )
  out <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(out)) {
      stop("Unknown argument: --", name, call. = FALSE)
    }
    if (i == length(argv)) {
      stop("Missing value for --", name, call. = FALSE)
    }
    value <- argv[[i + 1]]
    if (name %in% c("batch_size", "smooth_batch_size", "seed")) {
      value <- as.integer(value)
    }
    out[[name]] <- value
    i <- i + 2
  }
  out
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest)) {
    return(invisible(dest))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  invisible(dest)
}

ensure_cytotrace2_source <- function(source_cache) {
  r_dir <- file.path(source_cache, "R")
  extdata_dir <- file.path(source_cache, "extdata")
  dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)

  base_raw <- "https://raw.githubusercontent.com/digitalcytometry/cytotrace2/main/cytotrace2_r"
  files_r <- c("preprocess.R", "prediction.R", "postprocessing.R", "cytotrace2.R")
  files_ext <- c(
    "features_model_training_17.csv",
    "mt_dict_human_to_mouse.csv",
    "mt_human_alias.csv",
    "mt_mouse_alias.csv",
    "parameter_dict_19.rds"
  )

  for (name in files_r) {
    download_if_missing(
      sprintf("%s/R/%s", base_raw, name),
      file.path(r_dir, name)
    )
  }
  for (name in files_ext) {
    download_if_missing(
      sprintf("%s/inst/extdata/%s", base_raw, name),
      file.path(extdata_dir, name)
    )
  }
  list(r_dir = r_dir, extdata_dir = extdata_dir)
}

rewrite_source_text <- function(text, extdata_dir) {
  replacements <- c(
    "features_model_training_17.csv",
    "mt_dict_human_to_mouse.csv",
    "mt_human_alias.csv",
    "mt_mouse_alias.csv",
    "parameter_dict_19.rds"
  )
  for (fname in replacements) {
    pattern <- sprintf('system.file\\("extdata", "%s",\\s*package = "CytoTRACE2"\\)', gsub("\\.", "\\\\.", fname))
    text <- gsub(pattern, sprintf('file.path("%s", "%s")', extdata_dir, fname), text, perl = TRUE)
  }
  text
}

source_cytotrace2_functions <- function(source_cache) {
  paths <- ensure_cytotrace2_source(source_cache)
  files <- c("preprocess.R", "prediction.R", "postprocessing.R", "cytotrace2.R")
  for (name in files) {
    text <- paste(readLines(file.path(paths$r_dir, name), encoding = "UTF-8", warn = FALSE), collapse = "\n")
    text <- rewrite_source_text(text, normalizePath(paths$extdata_dir, winslash = "/", mustWork = TRUE))
    eval(parse(text = text), envir = .GlobalEnv)
  }
}

run_one <- function(row, species, batch_size, smooth_batch_size, seed) {
  message("RUN ", row$mode, " ", row$study_sample)
  scores <- cytotrace2(
    input = row$expression_tsv_gz,
    species = species,
    is_seurat = FALSE,
    batch_size = batch_size,
    smooth_batch_size = smooth_batch_size,
    parallelize_models = FALSE,
    parallelize_smoothing = FALSE,
    ncores = 1,
    seed = seed
  )
  scores$cell_id <- rownames(scores)
  cell_map <- utils::read.delim(row$cell_map_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  out <- merge(cell_map, scores, by = "cell_id", all.x = TRUE, sort = FALSE)
  out$mode <- row$mode
  out$study_sample <- row$study_sample
  out
}

params <- parse_args(args)
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
source_cytotrace2_functions(params$source_cache)

manifest <- utils::read.delim(params$input_manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!"mode" %in% names(manifest)) {
  stop("Input manifest must include a mode column.", call. = FALSE)
}
if (params$mode != "all") {
  manifest <- manifest[manifest$mode == params$mode, , drop = FALSE]
}
if (nrow(manifest) == 0) {
  stop("No CytoTRACE2 input rows selected.", call. = FALSE)
}

results <- lapply(split(manifest, seq_len(nrow(manifest))), function(part) {
  run_one(part[1, , drop = FALSE], params$species, params$batch_size, params$smooth_batch_size, params$seed)
})
combined <- do.call(rbind, results)

out_name <- if (params$mode == "all") "figure1c_cytotrace2_scores_by_cell.tsv.gz" else sprintf("figure1c_cytotrace2_scores_by_cell.%s.tsv.gz", params$mode)
out_path <- file.path(params$output_dir, out_name)
con <- gzfile(out_path, open = "wt")
utils::write.table(combined, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
message("WROTE ", normalizePath(out_path, winslash = "/", mustWork = FALSE))
