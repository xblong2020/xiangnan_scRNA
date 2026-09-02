#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_figure5_script_contract.R"
root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
scripts_dir <- file.path(root, "scripts")

required <- c(
  "figure5_plot_theme.R",
  "figure5_temporal_core.R",
  "figure5_00_preflight_audit.R",
  "figure5_00b_frozen_signature_audit.R",
  "figure5_01_calculate_three_axis_scores.R",
  "figure5_02_orient_pseudotime.R",
  "figure5_03_build_patient_pseudobulk.R",
  "plot_figure5a_temporal_framework.R",
  "figure5_04_fit_three_axis_gam.R",
  "plot_figure5b_three_axis_gam.R",
  "figure5_05_prepare_temporal_heatmap.R",
  "plot_figure5c_temporal_heatmap.R",
  "figure5_06_calculate_temporal_landmarks.R",
  "figure5_07_bootstrap_temporal_landmarks.R",
  "plot_figure5d_temporal_landmarks.R",
  "figure5_08_analyze_precedence_probability.R",
  "plot_figure5e_precedence_matrix.R",
  "figure5_09_analyze_method_concordance.R",
  "plot_figure5f_method_concordance.R",
  "figure5_10_analyze_patient_temporal_order.R",
  "plot_figure5g_patient_forest.R",
  "plot_figure5h_overlapping_phase_model.R",
  "plot_figure5_temporal_positioning_preview.R",
  "validate_figure5_temporal_positioning.R",
  "generate_figure5_temporal_positioning_report.R",
  "run_figure5_temporal_positioning.ps1"
)

missing <- required[!file.exists(file.path(scripts_dir, required))]
if (length(missing)) stop("Missing Figure 5 scripts: ", paste(missing, collapse = ", "), call. = FALSE)

plot_scripts <- required[grepl("^plot_.*\\.R$", required)]
for (script in plot_scripts) {
  text <- paste(readLines(file.path(scripts_dir, script), warn = FALSE), collapse = "\n")
  if (!grepl("figure5_plot_theme\\.R", text)) stop(script, " does not source the shared theme", call. = FALSE)
  if (grepl("matplotlib|seaborn|scanpy\\.pl|viridis|rainbow|\\bjet\\b", text, ignore.case = TRUE)) {
    stop(script, " contains a forbidden plotting tool or palette", call. = FALSE)
  }
}

source(file.path(scripts_dir, "figure5_plot_theme.R"), local = FALSE)
expected <- c(
  identity_loss = lancet_palette[[1]],
  stress_transition = lancet_palette[[3]],
  sox4_stabilization = lancet_palette[[2]]
)
if (!identical(unname(axis_palette), unname(expected))) stop("Axis palette does not match the Lancet contract", call. = FALSE)
if (!all(c("Supported", "Partial", "Not resolved", "Opposite", "Not available") %in% names(evidence_palette))) {
  stop("Evidence palette is incomplete", call. = FALSE)
}
if (!is.function(export_figure5_plot)) stop("Unified Figure 5 export function is missing", call. = FALSE)

h_text <- paste(readLines(file.path(scripts_dir, "plot_figure5h_overlapping_phase_model.R"), warn = FALSE), collapse = "\n")
if (grepl("geom_rect", h_text, fixed = TRUE)) stop("Figure 5H still contains hard rectangular stage boundaries", call. = FALSE)
if (!grepl("geom_polygon", h_text, fixed = TRUE) || !grepl("scale_alpha_identity", h_text, fixed = TRUE)) {
  stop("Figure 5H activity bands are not rendered as alpha-gradient polygons", call. = FALSE)
}
if (!grepl("uses_old_onset_for_formal_band = FALSE", h_text, fixed = TRUE) ||
    !grepl("onset_finite_fraction = onset", h_text, fixed = TRUE)) {
  stop("Figure 5H does not explicitly exclude old onset or audit corrected-onset stability", call. = FALSE)
}

cat("Figure 5 script contract tests passed\n")
