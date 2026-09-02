#!/usr/bin/env Rscript

## Shared data contracts and statistical helpers for Figure 6.

source(file.path(
  normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = FALSE),
  "scripts", "figure6_plot_theme.R"
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

FIGURE6_METADATA_DIR <- file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "figure6_directional_network")
FIGURE6_PROCESSED_DIR <- file.path(FIGURE6_PROJECT_ROOT, "data", "processed", "driver", "figure6_directional_network")
FIGURE6_REPORT_PATH <- file.path(FIGURE6_PROJECT_ROOT, "reports", "figure6_directional_network_report.md")

figure6_init_dirs <- function() {
  dirs <- c(
    FIGURE6_METADATA_DIR,
    FIGURE6_PROCESSED_DIR,
    file.path(FIGURE6_PROJECT_ROOT, "reports"),
    file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6_directional_network_preview")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

figure6_get_arg <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

figure6_norm_path <- function(path) {
  enc2utf8(gsub("\\\\", "/", normalizePath(path, mustWork = FALSE)))
}

figure6_write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    value, path, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null", digits = 16
  )
}

figure6_fread <- function(path, ...) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  data.table::fread(path, ...)
}

figure6_fwrite <- function(x, path, compress = grepl("\\.gz$", path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA", compress = if (compress) "gzip" else "none")
  invisible(path)
}

FIGURE6_AXIS_TFS <- list(
  identity_axis = c("HNF4A", "PPARA"),
  stress_axis = c("EGR1", "CEBPB", "JUN", "JUNB", "JUND", "FOS", "ATF3"),
  sox4_axis = "SOX4",
  control = c("HLF", "IRF1", "MAFB", "MAFF", "MYC")
)

FIGURE6_AP1_MEMBERS <- c("JUN", "JUNB", "JUND", "FOS", "ATF3")
FIGURE6_CORE_OUTPUTS <- c(
  "identity_program_change",
  "stress_transition_change",
  "sox4_programme_change",
  "malignant_fate_change",
  "proliferation_change",
  "cnv_malignant_signature_change"
)

FIGURE6_OUTPUT_LABELS <- c(
  identity_program_change = "Hepatocyte identity",
  stress_transition_change = "Stress-transition",
  sox4_programme_change = "SOX4 programme",
  malignant_fate_change = "Malignant fate",
  proliferation_change = "Proliferation",
  cnv_malignant_signature_change = "CNV-associated malignant\nexpression signature"
)

figure6_axis_for_tf <- function(tf) {
  ifelse(
    tf %in% FIGURE6_AXIS_TFS$identity_axis, "identity_axis",
    ifelse(
      tf %in% FIGURE6_AXIS_TFS$stress_axis, "stress_axis",
      ifelse(tf %in% FIGURE6_AXIS_TFS$sox4_axis, "sox4_axis", "control")
    )
  )
}

figure6_dataset_balanced_mean <- function(values, datasets) {
  ok <- is.finite(values) & !is.na(datasets)
  if (!any(ok)) return(NA_real_)
  means <- tapply(values[ok], as.character(datasets[ok]), mean, na.rm = TRUE)
  mean(means[is.finite(means)], na.rm = TRUE)
}

figure6_stratified_bootstrap <- function(sample_table, value_col, n_boot = 1000L, seed = 20260805L) {
  if (!nrow(sample_table) || !all(c(value_col, "dataset") %in% names(sample_table))) {
    return(list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, boot = numeric()))
  }
  work <- sample_table[is.finite(get(value_col)) & !is.na(dataset)]
  if (nrow(work) < 3L || uniqueN(work$dataset) < 1L) {
    return(list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, boot = numeric()))
  }
  estimate <- figure6_dataset_balanced_mean(work[[value_col]], work$dataset)
  set.seed(seed)
  by_dataset <- split(seq_len(nrow(work)), as.character(work$dataset))
  boot_values <- replicate(n_boot, {
    ix <- unlist(lapply(by_dataset, function(z) sample(z, length(z), replace = TRUE)), use.names = FALSE)
    figure6_dataset_balanced_mean(work[[value_col]][ix], work$dataset[ix])
  })
  boot_values <- boot_values[is.finite(boot_values)]
  if (!length(boot_values)) {
    return(list(estimate = estimate, ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, boot = numeric()))
  }
  pvalue <- 2 * min(mean(boot_values <= 0), mean(boot_values >= 0))
  if (pvalue == 0) pvalue <- 1 / length(boot_values)
  list(
    estimate = estimate,
    ci_low = unname(stats::quantile(boot_values, 0.025, na.rm = TRUE)),
    ci_high = unname(stats::quantile(boot_values, 0.975, na.rm = TRUE)),
    pvalue = min(pvalue, 1),
    boot = boot_values
  )
}

figure6_significance <- function(q) {
  ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
}

figure6_package_versions <- function(packages) {
  data.table(
    package = packages,
    available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(
      packages,
      function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_,
      character(1)
    )
  )
}

figure6_safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

figure6_jaccard <- function(a, b) {
  a <- unique(na.omit(as.character(a)))
  b <- unique(na.omit(as.character(b)))
  den <- length(union(a, b))
  if (!den) return(NA_real_)
  length(intersect(a, b)) / den
}

figure6_init_dirs()
