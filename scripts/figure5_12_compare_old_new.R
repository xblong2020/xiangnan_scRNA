#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
old_paths <- figure5_paths(root, create = FALSE)
new_paths <- figure5_onset_fix_paths(root)

summarize_bootstrap <- function(bootstrap) {
  bootstrap[, {
    finite <- time[is.finite(time)]
    list(
      bootstrap_median = if (length(finite)) stats::median(finite) else NA_real_,
      bootstrap_q025 = if (length(finite)) unname(stats::quantile(finite, 0.025, type = 8)) else NA_real_,
      bootstrap_q975 = if (length(finite)) unname(stats::quantile(finite, 0.975, type = 8)) else NA_real_,
      finite_fraction = length(finite) / .N
    )
  }, by = .(axis, landmark)]
}

old_landmarks <- fread(file.path(old_paths$metadata, "figure5_temporal_landmarks.tsv"))[, .(axis, landmark, old_estimate = time)]
new_landmarks <- fread(file.path(new_paths$metadata, "figure5_temporal_landmarks.tsv"))[, .(axis, landmark, new_estimate = time)]
old_boot <- as.data.table(readRDS(file.path(old_paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
new_boot <- as.data.table(readRDS(file.path(new_paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
old_summary <- summarize_bootstrap(old_boot)
setnames(old_summary, c("bootstrap_median", "bootstrap_q025", "bootstrap_q975", "finite_fraction"),
         c("old_bootstrap_median", "old_q025", "old_q975", "old_finite_fraction"))
new_summary <- summarize_bootstrap(new_boot)
setnames(new_summary, c("bootstrap_median", "bootstrap_q025", "bootstrap_q975", "finite_fraction"),
         c("new_bootstrap_median", "new_q025", "new_q975", "new_finite_fraction"))

comparison <- Reduce(function(x, y) merge(x, y, by = c("axis", "landmark"), all = TRUE),
                     list(old_landmarks, new_landmarks, old_summary, new_summary))
comparison[, absolute_change := abs(new_estimate - old_estimate)]
comparison[, old_ci := ifelse(is.finite(old_q025) & is.finite(old_q975), sprintf("[%.3f, %.3f]", old_q025, old_q975), "NA")]
comparison[, new_ci := ifelse(is.finite(new_q025) & is.finite(new_q975), sprintf("[%.3f, %.3f]", new_q025, new_q975), "NA")]
comparison[, interpretation_changed := fifelse(
  xor(is.finite(old_estimate), is.finite(new_estimate)) |
    xor(old_finite_fraction >= 0.80, new_finite_fraction >= 0.80) |
    (is.finite(absolute_change) & absolute_change > figure5_temporal_parameters$precedence_tolerance),
  "yes", "no"
)]
setcolorder(comparison, c(
  "axis", "landmark", "old_estimate", "new_estimate", "absolute_change",
  "old_bootstrap_median", "new_bootstrap_median", "old_ci", "new_ci",
  "old_finite_fraction", "new_finite_fraction", "interpretation_changed",
  "old_q025", "old_q975", "new_q025", "new_q975"
))
figure5_write_tsv(comparison, file.path(new_paths$metadata, "figure5_landmarks_old_vs_new.tsv"))

old_precedence <- fread(file.path(old_paths$metadata, "figure5e_precedence_probabilities.tsv"))
new_precedence <- fread(file.path(new_paths$metadata, "figure5e_precedence_probabilities.tsv"))
old_precedence <- old_precedence[, .(
  landmark, comparison, old_probability = probability,
  old_interval = ifelse(is.finite(ci_lower) & is.finite(ci_upper), sprintf("[%.3f, %.3f]", ci_lower, ci_upper), "NA"),
  old_median_delta = median_time_delta,
  old_evidence_grade = evidence_grade,
  old_effective_n = n_bootstrap_effective
)]
new_precedence <- new_precedence[, .(
  landmark, comparison, new_probability = probability,
  new_delta_interval = ifelse(is.finite(delta_q025) & is.finite(delta_q975), sprintf("[%.3f, %.3f]", delta_q025, delta_q975), "NA"),
  new_median_delta = median_delta,
  new_tie_fraction = tie_fraction,
  new_evidence_grade = evidence_grade,
  new_valid_n = n_valid,
  new_Monte_Carlo_SE = Monte_Carlo_SE
)]
precedence_comparison <- merge(old_precedence, new_precedence, by = c("landmark", "comparison"), all = TRUE)
precedence_comparison[, probability_change := new_probability - old_probability]
precedence_comparison[, interpretation_changed := fifelse(
  xor(is.finite(old_probability), is.finite(new_probability)) |
    (is.finite(probability_change) & abs(probability_change) > figure5_temporal_parameters$precedence_tolerance) |
    old_evidence_grade != new_evidence_grade,
  "yes", "no"
)]
precedence_comparison[, `:=`(
  landmark_order = match(landmark, c("onset", "maximum_slope", "t50")),
  comparison_order = match(comparison, precedence_pairs$comparison)
)]
setorder(precedence_comparison, landmark_order, comparison_order)
precedence_comparison[, c("landmark_order", "comparison_order") := NULL]
figure5_write_tsv(precedence_comparison, file.path(new_paths$metadata, "figure5_precedence_old_vs_new.tsv"))

original_onset_sensitivity <- comparison[landmark == "onset_time", .(
  axis, old_estimate, old_bootstrap_median, old_ci, old_finite_fraction,
  status = "sensitivity/source-data only; excluded from corrected Figure 5 main panels",
  old_definition = "first derivative >=10% of global maximum positive derivative, searched from pseudotime 0"
)]
figure5_write_tsv(original_onset_sensitivity, file.path(new_paths$metadata, "figure5_original_onset_sensitivity.tsv"))

old_predictions <- as.data.table(readRDS(file.path(old_paths$processed, "figure5b_gam_predictions.rds")))[scenario == "primary",
  .(axis, pseudotime, old_fit = fit)]
new_predictions <- as.data.table(readRDS(file.path(new_paths$processed, "figure5b_gam_predictions.rds")))[scenario == "primary",
  .(axis, pseudotime, new_fit = fit)]
curve_comparison <- merge(old_predictions, new_predictions, by = c("axis", "pseudotime"))
curve_summary <- curve_comparison[, .(
  max_absolute_difference = max(abs(new_fit - old_fit), na.rm = TRUE),
  rmse = sqrt(mean((new_fit - old_fit)^2, na.rm = TRUE)),
  pearson_correlation = stats::cor(old_fit, new_fit, use = "complete.obs"),
  old_range = diff(range(old_fit, na.rm = TRUE)),
  new_range = diff(range(new_fit, na.rm = TRUE))
), by = axis]
curve_summary[, reason := "coverage correction changed the primary inference unit from incomplete eligible patient tokens to all sample tokens"]
figure5_write_tsv(curve_summary, file.path(new_paths$metadata, "figure5b_gam_old_vs_new_curve.tsv"))

old_c <- fread(file.path(old_paths$metadata, "figure5c_temporal_landmarks.tsv"))[, .(entity, axis, old_maximum_slope_time = maximum_slope_time, old_display_order = display_order)]
new_c <- fread(file.path(new_paths$metadata, "figure5c_temporal_landmarks.tsv"))[, .(entity, axis, new_maximum_slope_time = maximum_slope_time, new_display_order = display_order)]
c_comparison <- merge(old_c, new_c, by = c("entity", "axis"), all = TRUE)
c_comparison[, `:=`(
  slope_absolute_change = abs(new_maximum_slope_time - old_maximum_slope_time),
  display_order_change = new_display_order - old_display_order
)]
figure5_write_tsv(c_comparison, file.path(new_paths$metadata, "figure5c_order_old_vs_new.tsv"))

figure5_write_json(list(
  landmark_rows = nrow(comparison),
  landmark_interpretations_changed = comparison[interpretation_changed == "yes", .N],
  precedence_rows = nrow(precedence_comparison),
  precedence_interpretations_changed = precedence_comparison[interpretation_changed == "yes", .N],
  original_onset_retained = TRUE,
  original_results_overwritten = FALSE
), file.path(new_paths$metadata, "figure5_old_vs_new_report.json"))
