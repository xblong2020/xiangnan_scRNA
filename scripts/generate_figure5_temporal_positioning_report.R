#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
paths <- figure5_paths(root)

preflight <- jsonlite::read_json(file.path(paths$metadata, "figure5_preflight_report.json"), simplifyVector = TRUE)
score_report <- jsonlite::read_json(file.path(paths$metadata, "figure5_three_axis_scores_report.json"), simplifyVector = TRUE)
orientation <- fread(file.path(paths$metadata, "figure5_pseudotime_orientation_audit.tsv"))
landmarks <- fread(file.path(paths$metadata, "figure5d_temporal_landmarks_plot_data.tsv"))
precedence <- fread(file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))[is_primary == TRUE]
concordance <- fread(file.path(paths$metadata, "figure5f_method_concordance.tsv"))
meta <- fread(file.path(paths$metadata, "figure5g_patient_meta_analysis.tsv"))[stratum == "overall"]
validation <- fread(file.path(paths$metadata, "figure5_validation_report.tsv"))
audit <- fread(file.path(paths$metadata, "figure5_frozen_signature_audit.tsv"))
hsummary <- fread(file.path(paths$metadata, "figure5h_activity_band_summary.tsv"))
hreport <- jsonlite::read_json(file.path(paths$metadata, "figure5h_overlapping_phase_model_report.json"), simplifyVector = TRUE)
conclusion <- figure5_result_text(precedence)
strongest <- conclusion
conservative <- "The three regulatory programmes occupied partially overlapping regions of the hepatocyte-state continuum. Stress-transition activity tended to precede SOX4-associated malignant-state stabilization, whereas the relative positioning of hepatocyte identity loss was less consistently resolved across methods."

fmt_table <- function(dt) {
  if (!nrow(dt)) return("Not available")
  header <- paste(names(dt), collapse = " | ")
  sep <- paste(rep("---", ncol(dt)), collapse = " | ")
  rows <- apply(dt, 1, function(x) paste(x, collapse = " | "))
  paste(c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |")), collapse = "\n")
}

package_versions <- figure5_package_versions(c("ggplot2", "ggsci", "patchwork", "scales", "data.table", "dplyr", "mgcv", "metafor", "jsonlite", "UCell", "reticulate"))
axis_counts <- audit[, .(n_genes = uniqueN(gene)), by = axis]
review_risks <- c(
  "No frozen patient_id field exists; explicit sample-name tokens were used and aggregate objects were excluded from patient meta-analysis.",
  "DPT and CellRank pseudotime were unavailable and remain labelled Not available.",
  "Cross-sectional pseudotime establishes relative positioning rather than physical time or causality.",
  "Figure 5H formal band boundaries exclude boundary-sensitive onset_time; directional activation fallback intervals are broader when t10 precedes supported high-score activation.",
  if (any(validation$status == "WARN")) paste("Validation warnings:", paste(validation[status == "WARN"]$check, collapse = "; ")) else NULL
)

report <- c(
  "# Figure 5 Temporal Positioning Report", "",
  "## Analysis identity", "",
  paste0("- R version: ", R.version.string),
  paste0("- Package versions: ", paste(names(package_versions), unlist(package_versions), sep = "=", collapse = "; ")),
  paste0("- Lancet palette HEX: ", paste(lancet_palette, collapse = ", ")),
  paste0("- Axis colours: identity loss=", axis_palette[[1]], "; stress transition=", axis_palette[[2]], "; SOX4 stabilization=", axis_palette[[3]]),
  paste0("- Random seed: 20260805"), "",
  "## Data objects and frozen signatures", "",
  paste0("Frozen strict cells: ", preflight$n_strict_cells, "; datasets: ", preflight$n_datasets, "; study samples: ", preflight$n_study_samples,
         "; CNV sample labels: ", preflight$n_cnv_samples, "; patient tokens: ", preflight$n_patient_tokens, "."),
  paste0("Three-axis signature gene counts: ", paste(axis_counts$axis, axis_counts$n_genes, sep = "=", collapse = "; "), "."),
  paste0("Primary scores: ", score_report$method, ". Standardization: ", score_report$standardization, "."), "",
  "## Pseudotime orientation audit", "", fmt_table(orientation[, .(method, status, flipped, oriented_score, start_state, end_state)]), "",
  "## Figure 5B GAM trends", "",
  "Patient-pseudotime-bin pseudobulk was fitted with mgcv GAMs using fixed dataset effects and a patient random-effect smooth where estimable. Confidence bands are 95% model-based intervals; cells are used for visualization and bin construction only.", "",
  "## Figure 5C temporal modules", "",
  "TF expression, cisTarget regulon AUC, frozen CellOracle/scTenifold targets and module scores were smoothed along oriented pseudotime, row-standardized, and ordered by maximum-slope time. Gene inclusion was fixed before temporal modelling.", "",
  "## Figure 5D temporal landmarks", "", fmt_table(landmarks[, .(axis, landmark, median, lower, upper, n_bootstrap)]), "",
  "## Figure 5E precedence probabilities", "", fmt_table(precedence[, .(comparison, probability, ci_lower, ci_upper, n_bootstrap_effective, evidence_grade)]), "",
  "## Figure 5F cross-method consistency", "", fmt_table(concordance[, .(method, comparison, probability, n_units, status)]), "",
  "## Figure 5G patient-level meta-analysis", "", fmt_table(meta[, .(comparison, status, n_patients, pooled_delta, ci_lower, ci_upper, I2, tau2, same_direction_fraction)]), "",
  "## Figure 5H overlapping-phase model", "",
  "Formal H band starts use bootstrap t10 when stable and directionally consistent; when t10 precedes sustained high-score activation, the first derivative-supported activation point is used. The original onset_time remains in source-data and sensitivity outputs but is excluded from the formal H band start.",
  fmt_table(hsummary[, .(axis, onset_bootstrap_median, t10, t10_lower, t10_upper, t50, maximum_slope, directional_activation, boundary_start, boundary_start_lower, boundary_start_upper, boundary_method, boundary_status, fade_start, fade_end)]),
  "Uncertainty is shown with light dashed outlines. Right edges are prominence-weighted fades and do not estimate a discrete programme end; background states are shown as observed median markers rather than hard stage rectangles.",
  "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events.", "",
  "## Sensitivity analyses", "",
  "Completed analyses include patient-pseudobulk, sample-balanced, dataset-balanced, CNV-strict, high-proliferation exclusion, proliferation adjustment, generic-stress removal, regulon AUC, TF expression, leave-one-dataset-out, leave-one-sample-out, available trajectory methods, k=4/5/6, and onset/t50/maximum-slope landmarks.", "",
  "## Unsupported or unresolved relationships", "",
  paste0("Primary conclusion selected by the prespecified thresholds: ", conclusion), "",
  "## Review-risk flags", "", paste0("- ", review_risks), "",
  "## SCI main-figure assessment", "",
  paste0("Validation status: ", if (any(validation$status == "FAIL")) "failed" else if (any(validation$status == "WARN")) "technically complete with review risks" else "technical requirements passed", "."), "",
  "## Recommended Results subheading", "",
  "Pseudotemporal positioning reveals overlapping phases of hepatocyte identity loss, stress-transition activation and SOX4-associated malignant-state stabilization.", "",
  "## Recommended Results paragraph", "",
  paste0(conclusion, " Patient-stratified bootstrap, available trajectory methods and random-effects patient analysis were used to quantify uncertainty. The programmes showed partially overlapping active intervals, so the result is interpreted as model-supported relative positioning along a cross-sectional hepatocyte-state continuum."), "",
  "## Conservative conclusion", "", conservative, "",
  "## Strongest data-supported conclusion", "", strongest, "",
  "## Claims that are not supported", "",
  "The analysis does not support real-time evolution, a definitive causal sequence, an obligatory sequential cascade, or experimentally established temporal order.", "",
  "## Recommended Figure legend", "",
  "Figure 5 | Temporal positioning reveals overlapping phases of hepatocyte identity loss, stress activation and malignant-state stabilization. (A) Analysis framework integrating oriented pseudotime, CNV/CellRank evidence and patient-level resampling. (B) Patient-pseudobulk GAM trends with 95% confidence intervals. (C) Temporally ordered frozen TF, regulon and target-gene heatmap. (D) Bootstrap median temporal landmarks and 95% intervals. (E) Precedence probabilities from stratified bootstrap. (F) Concordance across available trajectory methods and sensitivity analyses; unavailable methods are shown explicitly. (G) Patient-level maximum-slope differences and random-effects pooled estimates. (H) Smoothed-score activity bands with bootstrap-informed directional boundaries; the bands represent relative programme prominence and do not indicate discrete activation or termination events. Pseudotemporal positioning does not establish physical time or direct causality."
)
writeLines(report, file.path(paths$reports, "figure5_temporal_positioning_report.md"), useBytes = TRUE)

total <- list(
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE), R_version = R.version.string, packages = package_versions,
  lancet_palette = lancet_palette, axis_palette = axis_palette, signature_gene_counts = as.list(setNames(axis_counts$n_genes, axis_counts$axis)),
  counts = preflight[c("n_strict_cells", "n_datasets", "n_study_samples", "n_cnv_samples", "n_patient_tokens")],
  temporal_landmarks = split(landmarks, landmarks$axis), precedence = split(precedence, precedence$comparison),
  patient_meta_analysis = split(meta, meta$comparison), conclusion = conclusion, conservative_conclusion = conservative,
  figure5h = list(summary = hsummary, report = hreport),
  review_risk_flags = review_risks, validation_failures = validation[status == "FAIL"]$check,
  figure1_4_unchanged = validation[check_id == 25]$status == "PASS"
)
figure5_write_json(total, file.path(paths$metadata, "figure5_total_report.json"))

cat("Figure 5 report:", file.path(paths$reports, "figure5_temporal_positioning_report.md"), "\n")
cat("R version:", R.version.string, "\n")
cat("ggsci version:", as.character(packageVersion("ggsci")), "\n")
cat("Lancet palette:", paste(lancet_palette, collapse = ","), "\n")
cat("Counts:", paste(names(total$counts), unlist(total$counts), sep = "=", collapse = "; "), "\n")
cat("Precedence:", paste(precedence$comparison, sprintf("%.3f", precedence$probability), sep = "=", collapse = "; "), "\n")
cat("Conclusion:", conclusion, "\n")
cat("Figure 1-4 unchanged:", total$figure1_4_unchanged, "\n")
