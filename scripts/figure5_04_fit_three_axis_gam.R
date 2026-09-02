#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
set.seed(20260805)
pseudo <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
selection <- select_figure5_primary_pseudobulk(pseudo)
primary <- selection$data
analysis_unit_note <- selection$analysis_unit_note
if (nrow(primary) < 8L) stop("Insufficient patient-level pseudobulk rows for GAM")

scenarios <- list(
  primary = list(data = primary, suffix = "", k = 5L, adjusted = FALSE),
  adjusted_proliferation = list(data = primary, suffix = "", k = 5L, adjusted = TRUE),
  no_high_proliferation = list(data = primary[!is.finite(proliferation_score) | proliferation_score <= quantile(proliferation_score, 0.90, na.rm = TRUE)], suffix = "", k = 5L, adjusted = FALSE),
  no_generic_stress = list(data = primary, suffix = "no_generic", k = 5L, adjusted = FALSE),
  regulon_auc = list(data = primary, suffix = "regulon_auc", k = 5L, adjusted = FALSE),
  tf_expression = list(data = primary, suffix = "tf_expression", k = 5L, adjusted = FALSE),
  no_cell_cycle = list(data = primary, suffix = "no_cell_cycle", k = 5L, adjusted = FALSE),
  cnv_strict = list(data = primary[cnv_strict_fraction >= 0.5 | pseudotime <= 0.2], suffix = "", k = 5L, adjusted = FALSE),
  k4 = list(data = primary, suffix = "", k = 4L, adjusted = FALSE),
  k6 = list(data = primary, suffix = "", k = 6L, adjusted = FALSE),
  sample_balanced = list(data = pseudo[aggregation_unit == "sample" & method == "main/consensus pseudotime" & patient_meta_eligible == TRUE], suffix = "", k = 5L, adjusted = FALSE),
  dataset_balanced = list(data = primary[, .SD[sample(.N, min(.N, min(primary[, .N, by = dataset_id]$N)))], by = dataset_id], suffix = "", k = 5L, adjusted = FALSE)
)

predictions <- list()
summaries <- list()
models <- list()
for (scenario in names(scenarios)) {
  cfg <- scenarios[[scenario]]
  for (axis in names(figure5_axis_score_columns)) {
    base_col <- figure5_axis_score_columns[[axis]]
    score_col <- if (nzchar(cfg$suffix) && !(cfg$suffix == "no_generic" && axis != "stress_transition")) paste0(base_col, "_", cfg$suffix) else base_col
    if (!score_col %in% names(cfg$data) || sum(is.finite(cfg$data[[score_col]])) < 8L) {
      summaries[[paste(scenario, axis)]] <- data.table(scenario, axis, score_column = score_col, status = "Not available", n_rows = 0L)
      next
    }
    fit <- landmarks_from_table(cfg$data, score_col, axis, k = cfg$k, adjusted = cfg$adjusted, use_random_effect = TRUE)
    key <- paste(scenario, axis, sep = "::")
    models[[key]] <- fit$model
    pred <- copy(fit$predictions)
    pred[, `:=`(scenario = scenario, axis = axis, score_column = score_col, lower = fit - 1.96 * se, upper = fit + 1.96 * se)]
    predictions[[key]] <- pred
    if (is.null(fit$model)) {
      summaries[[key]] <- data.table(scenario, axis, score_column = score_col, status = "model_failed", n_rows = nrow(cfg$data))
    } else {
      sm <- summary(fit$model)
      srow <- if (nrow(sm$s.table)) sm$s.table[1, ] else rep(NA_real_, 4)
      summaries[[key]] <- data.table(scenario, axis, score_column = score_col, status = "fitted", n_rows = nrow(cfg$data),
                                      n_patients = uniqueN(cfg$data$patient_id), n_datasets = uniqueN(cfg$data$dataset_id),
                                      edf = unname(srow[[1]]), statistic = unname(srow[[3]]), p_value = unname(srow[[4]]),
                                      deviance_explained = sm$dev.expl, adjusted_r_squared = sm$r.sq,
                                      onset_time = fit$landmarks$onset_time, t10 = fit$landmarks$t10, t50 = fit$landmarks$t50,
                                      maximum_slope_time = fit$landmarks$maximum_slope_time, peak_time = fit$landmarks$peak_time,
                                      extremum_time = fit$landmarks$extremum_time, plateau_time = fit$landmarks$plateau_time,
                                      decline_onset = fit$landmarks$decline_onset,
                                      baseline_scale = fit$diagnostics$baseline_scale,
                                      baseline_scale_source = fit$diagnostics$baseline_scale_source,
                                      observed_coverage_min = fit$diagnostics$curve_coverage_min,
                                      observed_coverage_max = fit$diagnostics$curve_coverage_max,
                                      coverage_ok = fit$diagnostics$coverage_ok,
                                      onset_found = fit$diagnostics$onset_found,
                                      failure_reason = fit$diagnostics$failure_reason)
    }
  }
}

predictions <- rbindlist(predictions, fill = TRUE)
summaries <- rbindlist(summaries, fill = TRUE)
figure5_write_tsv(predictions, file.path(paths$metadata, "figure5b_gam_predictions.tsv.gz"))
figure5_write_tsv(summaries, file.path(paths$metadata, "figure5b_gam_model_summary.tsv"))
saveRDS(models, file.path(paths$processed, "figure5b_gam_models.rds"))
saveRDS(predictions, file.path(paths$processed, "figure5b_gam_predictions.rds"))
figure5_write_json(list(model = "mgcv::gam(score ~ s(pseudotime,k=5,bs='cr') + dataset_id + s(patient_id,bs='re'), method='REML')",
                        inference_unit = analysis_unit_note, random_seed = 20260805, scenarios = names(scenarios),
                        temporal_parameters = figure5_temporal_parameters,
                        baseline_scale = "raw patient/sample-pseudobulk baseline-window MAD, then baseline SD, then fixed 0.10 x fitted robust scale",
                        coverage_selection = list(selected_unit = selection$unit, observed_coverage = selection$coverage,
                                                  fallback_reason = selection$fallback_reason),
                        output_namespace = "figure5_temporal_positioning_onset_fix",
                        packages = figure5_package_versions(c("mgcv", "data.table", "jsonlite")), n_primary_rows = nrow(primary),
                        review_risk_flags = if (nzchar(selection$fallback_reason)) c(analysis_unit_note, selection$fallback_reason) else character()),
                   file.path(paths$metadata, "figure5b_gam_report.json"))
message("Figure 5B GAM fitting complete: ", nrow(summaries), " model rows.")
