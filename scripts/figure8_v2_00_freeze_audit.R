#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE))
  suppressWarnings(try(Sys.setlocale("LC_COLLATE", "Chinese"), silent = TRUE))
}

suppressPackageStartupMessages(library(data.table))

figure8_v2_project_root <- function() {
  explicit <- Sys.getenv("FIGURE8_PROJECT_ROOT", unset = "")
  if (nzchar(explicit)) return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(cwd, "AGENTS.md")) && dir.exists(file.path(cwd, "scripts"))) return(cwd)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
    return(normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

figure8_v2_relpath <- function(paths, root) {
  paths <- gsub("\\\\", "/", as.character(paths))
  root <- sub("/$", "", gsub("\\\\", "/", as.character(root)))
  prefix <- paste0(root, "/")
  ifelse(startsWith(paths, prefix), substring(paths, nchar(prefix) + 1L), paths)
}

figure8_v2_protected_files <- function(root = figure8_v2_project_root()) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  roots <- c("scripts", "metadata/driver", "figures/driver", "reports", "data/processed/driver")
  roots <- roots[dir.exists(roots)]
  files <- unlist(lapply(roots, list.files, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE), use.names = FALSE)
  files <- files[file.exists(files) & !dir.exists(files)]
  rel <- figure8_v2_relpath(files, root)

  excluded <- grepl(
    "(^|/)(figure8_transcriptomic_reversal_v2_mainfigure|figure8_v2_[^/]*|r_library|__pycache__|[.]pytest_cache)(/|$)",
    rel, ignore.case = TRUE, perl = TRUE
  ) |
    grepl("^scripts/(figure8_v2_|run_figure8_v2)", rel, ignore.case = TRUE) |
    grepl("^reports/figure8_transcriptomic_reversal_v2_mainfigure(/|$)", rel, ignore.case = TRUE)

  files <- files[!excluded]
  rel <- rel[!excluded]

  in_data <- startsWith(rel, "data/processed/driver/")
  protected_data <- grepl(
    "figure[1-7]|trajectory|celloracle|scenic|sctenifold|copykat|module[1-9]|driver_.*[.]h5ad|figure8_transcriptomic_reversal",
    rel, ignore.case = TRUE, perl = TRUE
  )
  keep <- !in_data | protected_data
  sort(unique(files[keep]))
}

figure8_v2_hash_manifest <- function(paths, root = figure8_v2_project_root(), workers = 4L) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  paths <- sort(unique(figure8_v2_relpath(paths, root)))
  if (!all(file.exists(paths))) stop("Cannot hash missing protected paths")
  workers <- max(1L, min(as.integer(workers), length(paths)))
  hash_one <- function(path) unname(tools::md5sum(path)[[1]])
  if (length(paths) && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    hashes <- unlist(parallel::parLapply(cl, paths, hash_one), use.names = FALSE)
  } else {
    hashes <- vapply(paths, hash_one, character(1))
  }
  info <- file.info(paths)
  data.table(
    file_path = figure8_v2_relpath(paths, root),
    size_bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    md5 = hashes
  )
}

figure8_v2_is_v1 <- function(rel) {
  grepl("^metadata/driver/figure8_transcriptomic_reversal(/|$)", rel, ignore.case = TRUE) |
    grepl("^data/processed/driver/figure8_transcriptomic_reversal(/|$)", rel, ignore.case = TRUE) |
    grepl("^figures/driver/figure8", rel, ignore.case = TRUE) |
    grepl("^reports/figure8_transcriptomic_reversal_report[.]md$", rel, ignore.case = TRUE) |
    grepl("^scripts/(figure8_|plot_figure8|validate_figure8|run_figure8)", rel, ignore.case = TRUE)
}

figure8_v2_design_changes <- function() {
  data.table(
    component = c(
      "DrugReflector input", "landmark QC", "cross-method validation", "orthogonal phenotype",
      "MoA/targets", "liver context", "matched-random specificity", "candidate ranking",
      "historical negative result", "main-figure decision"
    ),
    v1 = c(
      "sparse 150 UP + 150 DOWN primary",
      "47/300 landmark genes (15.7%)",
      "DrugReflector/L1000FWD/CLUE exact intersection; three-way overlap = 0",
      "unavailable",
      "compound-level annotations unavailable",
      "HEPG2/HCC515/HA1E treated as liver-relevant connectivity contexts",
      "1000 random signatures; specificity gate = FALSE",
      "single integrated heuristic plus evidence tiers",
      "1041 discordant candidates; no definitive consensus compound",
      "Extended Data recommendation"
    ),
    v2 = c(
      "continuous 978-landmark three-axis rescue v-score primary; sparse signature retained as sensitivity",
      "direct 978-model-gene audit with usable/non-zero and axis-balance metrics",
      "cross-framework corroboration across related LINCS/CMap resources; zero overlap retained",
      "PRISM 23Q2 primary plus 19Q4 secondary cancer-cell viability when mapped",
      "curated MoA/targets separated from Figure 6 network-consistent inference",
      "metadata-verified adult HCC, other liver cancer, hepatoblastoma, and non-liver contexts",
      "2000 expression/variance/detection/axis/absolute-weight matched null signatures with two-sided empirical P",
      "visible evidence-component matrix, conservative and coverage-aware summaries, frozen gates",
      "all v1 negatives retained as historical baseline",
      "automatic MAIN_FIGURE_READY or EXTENDED_DATA_ONLY"
    ),
    reason = c(
      "match DrugReflector input geometry and reduce sparse-projection bias",
      "measure information over the actual model space without a 25% target",
      "avoid overstating non-independent platform evidence",
      "add orthogonal phenotype without interpreting it as safety",
      "add mechanistic annotation while preserving provenance boundaries",
      "correct HCC515/HA1E context and retain HepG2 caveat",
      "test signature specificity against a stringent empirical null",
      "make every evidence contribution and missing layer auditable",
      "prevent selective deletion of negative findings",
      "let frozen evidence gates, not the desired narrative, determine placement"
    )
  )
}

figure8_v2_write_baseline <- function(root = figure8_v2_project_root(), workers = 4L) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  out_dir <- "metadata/driver/figure8_transcriptomic_reversal_v2_mainfigure"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  before_f17 <- file.path(out_dir, "figure8_v2_protected_figure1_7_hash_before.tsv")
  before_v1 <- file.path(out_dir, "figure8_v2_protected_figure8_v1_hash_before.tsv")
  design_path <- file.path(out_dir, "figure8_v1_vs_v2_design_changes.tsv")
  if (file.exists(before_f17) || file.exists(before_v1)) {
    stop("Protection baseline already exists; refusing to overwrite it. Remove only after an explicit baseline-reset decision.")
  }

  paths <- figure8_v2_protected_files(root)
  manifest <- figure8_v2_hash_manifest(paths, root, workers)
  is_v1 <- figure8_v2_is_v1(manifest$file_path)
  fwrite(manifest[!is_v1], before_f17, sep = "\t", quote = FALSE, na = "NA")
  fwrite(manifest[is_v1], before_v1, sep = "\t", quote = FALSE, na = "NA")
  fwrite(figure8_v2_design_changes(), design_path, sep = "\t", quote = FALSE, na = "NA")

  summary <- data.table(
    snapshot_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    figure1_7_and_frozen_files = sum(!is_v1),
    figure8_v1_files = sum(is_v1),
    total_size_gb = sum(manifest$size_bytes) / 1024^3,
    workers = as.integer(workers)
  )
  fwrite(summary, file.path(out_dir, "figure8_v2_protection_baseline_summary.tsv"), sep = "\t", quote = FALSE)
  invisible(list(figure1_7 = before_f17, figure8_v1 = before_v1, design = design_path, summary = summary))
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  workers_arg <- grep("^--workers=", args, value = TRUE)
  workers <- if (length(workers_arg)) as.integer(sub("^--workers=", "", workers_arg[[1]])) else 4L
  result <- figure8_v2_write_baseline(workers = workers)
  print(result$summary)
}
