#!/usr/bin/env Rscript

# Focused contract test for the refactored Figure 5 namespace.  The full
# validation script writes the machine-readable report; this test makes the
# high-risk panel-remapping and boundary rules independently executable.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root, create = FALSE)
validation <- fread(file.path(paths$metadata, "figure5_six_panel_validation_report.tsv"))
if (!nrow(validation) || !all(validation$passed)) stop("Figure 5 six-panel validation report contains failures", call. = FALSE)

summary <- fread(file.path(paths$metadata, "figure5_six_panel_F_activity_band_summary.tsv"))
if (!all(summary$boundary_method[summary$t10_stable] == "bootstrap_t10")) stop("5F does not enforce t10-priority starts", call. = FALSE)
if (any(summary$boundary_start == summary$onset, na.rm = TRUE)) stop("5F formal start equals onset_time", call. = FALSE)

evidence <- fread(file.path(paths$metadata, "figure5_six_panel_E_precedence_probabilities.tsv"))
if (!all(c("onset", "t10", "t50", "maximum_slope") %chin% evidence$landmark)) stop("5E landmark contract failed", call. = FALSE)
if (any(evidence$landmark %chin% c("peak", "peak_time"))) stop("Peak leaked into 5E", call. = FALSE)

cat("Figure 5 six-panel focused contract test passed\n")

