#!/usr/bin/env Rscript

## Shared Figure 8 paths, colour contract, typography, export, and small utilities.

figure8_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else ""
  root <- normalizePath(
    if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
    winslash = "/", mustWork = FALSE
  )
  if (!dir.exists(file.path(root, "scripts"))) {
    root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
  root
}

FIGURE8_PROJECT_ROOT <- figure8_project_root()
FIGURE8_R_LIBRARY <- file.path(
  FIGURE8_PROJECT_ROOT, "data", "processed", "driver",
  "figure8_transcriptomic_reversal", "r_library"
)
if (dir.exists(FIGURE8_R_LIBRARY)) .libPaths(c(FIGURE8_R_LIBRARY, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(jsonlite)
  library(patchwork)
  library(ggrepel)
})

FIGURE8_METADATA_DIR <- file.path(
  FIGURE8_PROJECT_ROOT, "metadata", "driver", "figure8_transcriptomic_reversal"
)
FIGURE8_DATA_DIR <- file.path(
  FIGURE8_PROJECT_ROOT, "data", "processed", "driver", "figure8_transcriptomic_reversal"
)
FIGURE8_REPORT_DIR <- file.path(FIGURE8_PROJECT_ROOT, "reports")
FIGURE8_SEED <- 20260805L

lancet_palette <- ggsci::pal_lancet("lanonc")(9)

axis_palette <- c(
  identity_rescue = lancet_palette[1],
  stress_suppression = lancet_palette[3],
  sox4_suppression = lancet_palette[2]
)

method_palette <- c(
  DrugReflector = lancet_palette[5],
  L1000FWD = lancet_palette[4],
  CLUE = lancet_palette[1],
  external_signature = lancet_palette[3]
)

evidence_palette <- c(
  tier_A = lancet_palette[3],
  tier_B = lancet_palette[4],
  tier_C = lancet_palette[6],
  exploratory = lancet_palette[8],
  discordant = lancet_palette[2],
  unresolved = lancet_palette[7],
  unavailable = "#FFFFFF"
)

reversal_gradient <- c(
  low = lancet_palette[2],
  mid = "#F7F7F7",
  high = lancet_palette[1]
)

workflow_palette <- c(
  completed = lancet_palette[3],
  partial = lancet_palette[6],
  unavailable = lancet_palette[8],
  failed = lancet_palette[2],
  external_validation = lancet_palette[1]
)

figure8_theme <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(size = 10, face = "plain", hjust = 0, colour = lancet_palette[9]),
      plot.subtitle = element_text(size = 8, colour = lancet_palette[9], hjust = 0),
      plot.caption = element_text(size = 6.8, colour = lancet_palette[9], hjust = 0),
      plot.tag = element_text(size = 10, face = "bold", family = "sans", colour = lancet_palette[9]),
      axis.title = element_text(size = 8.5, colour = lancet_palette[9]),
      axis.text = element_text(size = 7.5, colour = lancet_palette[9]),
      axis.line = element_line(linewidth = 0.4, colour = lancet_palette[9]),
      axis.ticks = element_line(linewidth = 0.4, colour = lancet_palette[9]),
      legend.title = element_text(size = 7.2, colour = lancet_palette[9]),
      legend.text = element_text(size = 7, colour = lancet_palette[9]),
      legend.key.height = grid::unit(0.34, "cm"),
      strip.text = element_text(size = 8, face = "bold", colour = lancet_palette[9]),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 6)
    )
}

figure8_dirs <- function() {
  dirs <- c(FIGURE8_METADATA_DIR, FIGURE8_DATA_DIR, FIGURE8_REPORT_DIR)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

figure8_panel_dir <- function(panel) {
  key <- tolower(panel)
  mapping <- c(
    a = "figure8a_target_state",
    b = "figure8b_signature_composition",
    c = "figure8c_reversal_workflow",
    d = "figure8d_drugreflector_stability",
    e = "figure8e_cross_method_concordance",
    f = "figure8f_mechanism_classes",
    g = "figure8g_external_signature_validation",
    h = "figure8h_integrated_prioritization",
    preview = "figure8_transcriptomic_reversal_preview"
  )
  if (!key %in% names(mapping)) stop("Unknown Figure 8 panel: ", panel)
  file.path(FIGURE8_PROJECT_ROOT, "figures", "driver", unname(mapping[[key]]))
}

figure8_save <- function(plot, panel, stem, width, height, save_rds = TRUE) {
  out_dir <- figure8_panel_dir(panel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- file.path(out_dir, paste0(stem, c(".pdf", ".png", ".svg", ".tiff")))
  names(outputs) <- c("pdf", "png", "svg", "tiff")
  ggplot2::ggsave(outputs[["pdf"]], plot = plot, device = grDevices::cairo_pdf,
                  width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
  ggplot2::ggsave(outputs[["png"]], plot = plot, device = "png", dpi = 600,
                  width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
  grDevices::svg(outputs[["svg"]], width = width, height = height, family = "sans", bg = "white")
  print(plot)
  grDevices::dev.off()
  ggplot2::ggsave(outputs[["tiff"]], plot = plot, device = "tiff", dpi = 600,
                  width = width, height = height, units = "in", compression = "lzw",
                  bg = "white", limitsize = FALSE)
  if (isTRUE(save_rds)) saveRDS(plot, file.path(FIGURE8_DATA_DIR, paste0(stem, "_plot.rds")))
  invisible(outputs)
}

figure8_write_tsv <- function(x, name, compress = FALSE) {
  figure8_dirs()
  suffix <- if (compress && !grepl("\\.gz$", name)) ".gz" else ""
  path <- file.path(FIGURE8_METADATA_DIR, paste0(name, suffix))
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
  path
}

figure8_write_data_tsv <- function(x, name, compress = FALSE) {
  figure8_dirs()
  suffix <- if (compress && !grepl("\\.gz$", name)) ".gz" else ""
  path <- file.path(FIGURE8_DATA_DIR, paste0(name, suffix))
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
  path
}

figure8_write_json <- function(x, name) {
  figure8_dirs()
  path <- file.path(FIGURE8_METADATA_DIR, name)
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  path
}

figure8_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

figure8_norm01 <- function(x, higher_is_better = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  rng <- range(x[ok])
  if (diff(rng) == 0) out[ok] <- 0.5 else out[ok] <- (x[ok] - rng[[1]]) / diff(rng)
  if (!higher_is_better) out[ok] <- 1 - out[ok]
  out
}

figure8_rank_score <- function(rank, n_compounds = 9597L) {
  rank <- suppressWarnings(as.numeric(rank))
  pmax(0, pmin(1, 1 - (rank - 1) / max(n_compounds - 1, 1)))
}

figure8_safe_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

figure8_package_versions <- function(packages) {
  setNames(lapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }), packages)
}

figure8_palette_contract <- function() {
  list(
    r_version = R.version.string,
    ggsci_version = as.character(utils::packageVersion("ggsci")),
    source = "ggsci::pal_lancet('lanonc')(9)",
    lancet_palette = as.list(unname(lancet_palette)),
    axis_palette = as.list(axis_palette),
    method_palette = as.list(method_palette),
    evidence_palette = as.list(evidence_palette),
    reversal_gradient = as.list(reversal_gradient),
    reversal_semantics = list(high = "desired rescue direction", midpoint = "no evident reversal", low = "disease direction or discordance")
  )
}

figure8_require <- function(paths, label = "required input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop(label, " missing: ", paste(missing, collapse = ", "))
  invisible(paths)
}

figure8_fread <- function(path, ...) {
  data.table::fread(path, na.strings = c("", "NA", "NaN", "null"), ...)
}

figure8_axis_from_component <- function(component) {
  component <- as.character(component)
  dplyr::case_when(
    component %in% c("hnf4a_ppara_rescue", "mature_hepatocyte", "tier1_rescue") ~ "identity_rescue",
    component %in% c("ap1_stress_proliferation", "cebpb_egr1_malignant_target") ~ "stress_suppression",
    component %in% c("sox4_state_specific", "c_malignant_like_fate") ~ "sox4_suppression",
    TRUE ~ "unresolved"
  )
}

figure8_root_path <- function(...) file.path(FIGURE8_PROJECT_ROOT, ...)

