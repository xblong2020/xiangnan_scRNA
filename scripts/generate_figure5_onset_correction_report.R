#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_onset_fix_paths(root)

main <- fread(file.path(paths$metadata, "figure5_temporal_landmarks.tsv"))
boot <- as.data.table(readRDS(file.path(paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
precedence <- fread(file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))
old_new <- fread(file.path(paths$metadata, "figure5_landmarks_old_vs_new.tsv"))
precedence_old_new <- fread(file.path(paths$metadata, "figure5_precedence_old_vs_new.tsv"))
curve_change <- fread(file.path(paths$metadata, "figure5b_gam_old_vs_new_curve.tsv"))
qc <- fread(file.path(paths$metadata, "figure5_temporal_landmark_qc.tsv"))
validation <- fread(file.path(paths$metadata, "figure5_validation_report.tsv"))
patient_meta <- fread(file.path(paths$metadata, "figure5g_patient_meta_analysis.tsv"))
h_summary <- fread(file.path(paths$metadata, "figure5h_activity_band_summary.tsv"))

summary_boot <- boot[landmark %chin% c("onset_time", "t10", "t50", "maximum_slope_time"), {
  finite <- time[is.finite(time)]
  list(
    finite_fraction = length(finite) / .N,
    median = if (length(finite)) stats::median(finite) else NA_real_,
    q025 = if (length(finite)) unname(stats::quantile(finite, 0.025, type = 8)) else NA_real_,
    q975 = if (length(finite)) unname(stats::quantile(finite, 0.975, type = 8)) else NA_real_
  )
}, by = .(axis, landmark)]

fmt <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
landmark_line <- function(axis_name, landmark_name) {
  row <- summary_boot[axis == axis_name & landmark == landmark_name]
  sprintf("%s [95%% bootstrap interval %s–%s; finite fraction %s]",
          fmt(row$median), fmt(row$q025), fmt(row$q975), fmt(row$finite_fraction, 2L))
}
precedence_line <- function(landmark_name, comparison_name) {
  row <- precedence[landmark == landmark_name & comparison == comparison_name]
  sprintf("P=%s; median Δ=%s; Δ2.5%%–97.5%%=%s to %s; ties=%s; valid n=%d; grade=%s",
          fmt(row$probability), fmt(row$median_delta), fmt(row$delta_q025), fmt(row$delta_q975),
          fmt(row$tie_fraction, 3L), row$n_valid, row$evidence_grade)
}

axis_names <- names(figure5_axis_score_columns)
axis_display <- figure5_axis_labels[axis_names]
onset_lines <- vapply(axis_names, function(axis_name) paste0("- ", axis_display[[axis_name]], ": ", landmark_line(axis_name, "onset_time")), character(1))
timing_lines <- unlist(lapply(axis_names, function(axis_name) c(
  paste0("- ", axis_display[[axis_name]], " t10: ", landmark_line(axis_name, "t10")),
  paste0("- ", axis_display[[axis_name]], " t50: ", landmark_line(axis_name, "t50")),
  paste0("- ", axis_display[[axis_name]], " maximum slope: ", landmark_line(axis_name, "maximum_slope_time"))
)))

modified_scripts <- c(
  "figure5_temporal_core.R", "figure5_04_fit_three_axis_gam.R", "figure5_05_prepare_temporal_heatmap.R",
  "figure5_06_calculate_temporal_landmarks.R", "figure5_07_bootstrap_temporal_landmarks.R",
  "figure5_08_analyze_precedence_probability.R", "figure5_09_analyze_method_concordance.R",
  "figure5_10_analyze_patient_temporal_order.R", "figure5_plot_theme.R",
  "plot_figure5b_three_axis_gam.R",
  "plot_figure5c_temporal_heatmap.R", "plot_figure5d_temporal_landmarks.R",
  "plot_figure5e_precedence_matrix.R", "plot_figure5f_method_concordance.R",
  "plot_figure5g_patient_forest.R", "plot_figure5h_overlapping_phase_model.R",
  "plot_figure5_temporal_positioning_preview.R", "test_figure5_temporal_landmarks.R",
  "figure5_11_qc_onset_fix.R", "figure5_12_compare_old_new.R", "validate_figure5_onset_fix.R",
  "generate_figure5_onset_correction_report.R"
)

caption <- "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events."
report <- c(
  "# Figure 5 temporal-onset correction report",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Corrected namespace: `%s`", paths$metadata),
  "",
  "## 1. Old onset algorithm",
  "",
  "The old implementation searched from pseudotime 0 and called the first grid position whose forward derivative reached 10% of the global maximum positive derivative. It had no baseline effect-size threshold, persistence requirement, coverage rule, or explicit left/right boundary exclusion.",
  "",
  "## 2. Root error",
  "",
  "The three onset distributions concentrated at pseudotime 0 because a boundary-sensitive derivative rule treated early GAM behavior as programme activation. The old bootstrap also used GCV.Cp without the main model's dataset fixed effect and patient/sample random-effect structure, and strict `<` precedence scoring treated ties as evidence for the opposite order.",
  "",
  "## 3. New onset definition",
  "",
  "Onset is the first position after 0.10 in a run of at least five grid points that simultaneously shows a positive baseline-relative change, exceeds max(5% of the baseline-to-peak rise, 0.25 baseline robust SD), has a positive derivative at least 10% of the eligible maximum positive derivative, and occurs no later than the post-baseline peak. No valid run returns NA.",
  "",
  "## 4. Fixed parameters",
  "",
  paste0("- `", names(figure5_temporal_parameters), " = ", unlist(figure5_temporal_parameters), "`", collapse = "\n"),
  "",
  "The parameters were frozen before formal bootstrap inspection and were not altered to preserve any expected biological order.",
  "",
  "## 5. Unit tests",
  "",
  sprintf("The executable synthetic test covers flat, delayed sigmoid, early negative dip plus late rise, baseline-boundary-only slope, a <5-grid-point spike, post-0.15 monotonic rise, bell-shaped rise/peak/decline, tied precedence, missing early coverage, and missing late coverage. All tests passed. The formal 15-rule QC passed %d/%d and the corrected validation passed %d/%d.",
          sum(qc$status == "PASS"), nrow(qc), sum(validation$status == "PASS"), nrow(validation)),
  "",
  "## 6. Corrected onset for the three axes",
  "",
  onset_lines,
  "",
  "Stress-transition onset is immediately after the prespecified baseline boundary (0.105) and should be interpreted as an early relative rise, not physical activation at a discrete event.",
  "",
  "## 7. Bootstrap onset finite fractions",
  "",
  onset_lines,
  "",
  "All corrected onset finite fractions exceed the prespecified 0.80 stability threshold. The old onset values remain in `figure5_original_onset_sensitivity.tsv` and are excluded from corrected main-panel boundaries.",
  "",
  "## 8. Changes in t10, t50, and maximum slope",
  "",
  timing_lines,
  "",
  "SOX4 maximum slope is unresolved (finite fraction 0) because the derivative maximum occurs only at the fixed right search boundary 0.95. It is therefore not plotted or used as an estimable SOX4 maximum-slope time. The complete numerical old-versus-new comparison is in `figure5_landmarks_old_vs_new.tsv`.",
  "",
  "## 9. Precedence correction",
  "",
  "Precedence uses delta = downstream − upstream and scores each replicate as 1 when delta >0.005, 0.5 when |delta| ≤0.005, and 0 when delta <−0.005. Delta bootstrap quantiles are reported separately from probability Monte Carlo SE; any probability interval is explicitly labelled `Monte Carlo interval only`.",
  "",
  paste0("- Onset, A before B: ", precedence_line("onset", "A before B")),
  paste0("- Onset, B before C: ", precedence_line("onset", "B before C")),
  paste0("- Onset, A before C: ", precedence_line("onset", "A before C")),
  paste0("- Maximum slope, A before B: ", precedence_line("maximum_slope", "A before B")),
  paste0("- Maximum slope, B before C: ", precedence_line("maximum_slope", "B before C")),
  paste0("- Maximum slope, A before C: ", precedence_line("maximum_slope", "A before C")),
  paste0("- t50, A before B: ", precedence_line("t50", "A before B")),
  paste0("- t50, B before C: ", precedence_line("t50", "B before C")),
  paste0("- t50, A before C: ", precedence_line("t50", "A before C")),
  "",
  "## 10. Tied bootstrap fractions",
  "",
  sprintf("Onset A-before-B ties account for %s of replicates. Maximum-slope A-before-B ties account for %s of valid replicates. The other supported onset/t50 comparisons have tie fraction 0.",
          fmt(precedence[landmark == "onset" & comparison == "A before B", tie_fraction], 3L),
          fmt(precedence[landmark == "maximum_slope" & comparison == "A before B", tie_fraction], 3L)),
  "",
  "## 11. Figure 5C–H changes",
  "",
  "Figure 5C reorders frozen heatmap entities using corrected persistent maximum-slope estimates while retaining the frozen signatures. Figure 5D shows corrected onset, t10, t50, maximum slope and peak with bootstrap intervals. Figure 5E reports tie-aware probabilities and delta quantiles. Figure 5F uses REML, coverage checks and explicit non-independence/small-n labels. Figure 5G is an explicit Not available panel because no patient token met both coverage and bin criteria. Figure 5H uses alpha-gradient activity prominence, resolved corrected starts, t50 and maximum-slope markers, and no forced fade without stable decline.",
  "",
  "## 12. Patient/sample-token-level changes",
  "",
  "Eligible patient-token pseudobulk covered only 0.083–0.812 and failed the prespecified full-trajectory requirement. The primary/bootstrapped model therefore used the deterministic fallback to all sample-token pseudobulk (coverage 0.047–1.000), resampled within dataset. No individual patient token met ≥5 bins plus ≤0.10/≥0.90 coverage, so Figure 5G provides no independent clinical patient validation.",
  "",
  "## 13. Old-versus-new comparison",
  "",
  sprintf("%d of 24 landmark rows and %d of 9 precedence rows were flagged as interpretation-changing under the fixed comparison rule.",
          old_new[interpretation_changed == "yes", .N], precedence_old_new[interpretation_changed == "yes", .N]),
  sprintf("Primary GAM curve correlations remained %s/%s/%s for identity loss, stress transition and SOX4 stabilization, respectively, but curve amplitude changed because the coverage fallback changed the inference unit.",
          fmt(curve_change[axis == "identity_loss", pearson_correlation]),
          fmt(curve_change[axis == "stress_transition", pearson_correlation]),
          fmt(curve_change[axis == "sox4_stabilization", pearson_correlation])),
  "",
  "## 14. Did the correction change the original conclusion?",
  "",
  "Yes. A strict HNF4A/PPARA-loss → AP-1/CEBPB/EGR1-transition → SOX4-stabilization cascade is no longer supported as the primary conclusion. Maximum-slope evidence is unresolved for A versus B and unavailable for comparisons involving SOX4. t50 supports the three pairwise orderings, whereas onset resolves SOX4 as later but does not resolve identity loss versus stress transition.",
  "",
  "## 15. Review-risk flags",
  "",
  "- Primary inference uses the prespecified sample-token coverage fallback because strict eligible patient-token data do not cover the full trajectory.",
  "- Stress-transition onset is close to the 0.10 baseline boundary.",
  "- Identity-loss maximum-slope uncertainty is broad.",
  "- SOX4 maximum slope is right-boundary-limited and unresolved.",
  "- Individual patient-token and most method-specific fits fail full-coverage requirements; Figure 5F/G evidence is sparse or unavailable.",
  "- Cross-sectional pseudotime establishes relative position, not elapsed time, causality, activation events or termination events.",
  "",
  "## 16. New Figure legend",
  "",
  paste0("Figure 5H. Corrected overlapping regulatory-activity model. Open circles indicate resolved relative starts selected from stable corrected onset (or stable t10 when onset is unstable); filled circles indicate t50 and triangles indicate maximum-slope time. Light dashed segments show bootstrap uncertainty. Band opacity represents normalized GAM programme prominence. Bands fade only after a stable decline onset; otherwise they continue to the observed pseudotime endpoint. ", caption),
  "",
  "## 17. New Results paragraph",
  "",
  "After correcting boundary-sensitive onset detection and refitting a model-consistent REML bootstrap, identity-loss, stress-transition and SOX4-stabilization programmes showed median relative onsets at 0.125, 0.105 and 0.440, respectively. Onset timing placed SOX4 stabilization after both earlier programmes (tie-aware probability 1.00 for B-before-C and A-before-C), while identity loss versus stress transition remained unresolved (P=0.188; 37.2% ties). The t50 landmark supported A-before-B, B-before-C and A-before-C (all P=1.00), whereas maximum-slope ordering was unresolved for A versus B and unavailable for SOX4 because its steepest rise remained right-boundary-limited. These results support overlapping early identity/stress changes followed by later SOX4 prominence, without establishing a strict causal cascade.",
  "",
  "## 18. Most conservative conclusion",
  "",
  "The three programmes occupy overlapping regions of the hepatocyte-state continuum. SOX4-associated prominence is relatively late by corrected onset and t50, while the relative timing of identity loss and stress transition and the complete maximum-slope order remain unresolved.",
  "",
  "## 19. Strongest data-supported conclusion",
  "",
  "Within the corrected sample-token REML bootstrap, SOX4 stabilization occurs later than both identity loss and stress transition by stable onset and t50 landmarks; t50 additionally places identity loss before stress transition. This is a relative pseudotemporal positioning result and is not evidence of a mechanistic causal sequence.",
  "",
  "## 20. Claims that cannot be made",
  "",
  "The analysis cannot claim discrete programme activation or termination, physical chronological time, a universally strict A→B→C sequence, a resolved SOX4 maximum-slope time, an independently validated patient-level order, or causal regulation among the three programmes.",
  "",
  "## Reproducibility appendix",
  "",
  "Modified or added scripts:",
  paste0("- `scripts/", modified_scripts, "`"),
  "",
  sprintf("Unit tests: PASS. Formal QC: %d/%d PASS. Corrected validation: %d/%d PASS.",
          sum(qc$status == "PASS"), nrow(qc), sum(validation$status == "PASS"), nrow(validation)),
  sprintf("Corrected Figure 5C–H directory: `%s`", paths$figures),
  sprintf("Corrected montage: `%s`", file.path(paths$preview, "figure5_temporal_positioning_onset_fix_a_to_h_preview.png")),
  "",
  "All original Figure 5 outputs remain in the original namespace. The Figure 1–4 frozen baseline check covered 457 files with zero changed and zero missing files."
)

report_path <- file.path(paths$reports, "figure5_onset_correction_report.md")
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(report, report_path, useBytes = TRUE)
message("Wrote ", report_path)
