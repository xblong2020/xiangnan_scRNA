#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
suppressPackageStartupMessages(library(metafor))
pseudo <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
dt <- pseudo[aggregation_unit == "patient" & method == "main/consensus pseudotime" & eligible_patient == TRUE]
if (uniqueN(dt$patient_id) < 2L) stop("Insufficient eligible patient tokens for meta-analysis")

comparison_map <- list(
  "A to B" = c("identity_loss", "stress_transition"),
  "B to C" = c("stress_transition", "sox4_stabilization"),
  "A to C" = c("identity_loss", "sox4_stabilization")
)

patient_axis_time <- function(data, score_col, axis_name) {
  work <- data[is.finite(pseudotime) & is.finite(get(score_col))]
  coverage <- if (nrow(work)) range(work$pseudotime, na.rm = TRUE) else c(NA_real_, NA_real_)
  if (nrow(work) < 5L || uniqueN(work$pseudotime) < 5L || !all(is.finite(coverage)) ||
      coverage[[1L]] > figure5_temporal_parameters$baseline_end || coverage[[2L]] < 0.90) return(NA_real_)
  k <- min(4L, uniqueN(work$pseudotime) - 1L)
  fit <- landmarks_from_table(work, score_col, axis_name, k = k,
                              adjusted = FALSE, use_random_effect = FALSE)
  fit$landmarks$maximum_slope_time
}

patient_rows <- list()
for (patient in unique(dt$patient_id)) {
  subset <- dt[patient_id == patient]
  n_bins <- uniqueN(subset$pseudotime_bin)
  coverage <- range(subset$pseudotime, na.rm = TRUE)
  coverage_ok <- n_bins >= 5L && coverage[[1L]] <= figure5_temporal_parameters$baseline_end && coverage[[2L]] >= 0.90
  times <- if (coverage_ok) {
    vapply(names(figure5_axis_score_columns), function(axis_name) {
      patient_axis_time(subset, figure5_axis_score_columns[[axis_name]], axis_name)
    }, numeric(1))
  } else {
    stats::setNames(rep(NA_real_, length(figure5_axis_score_columns)), names(figure5_axis_score_columns))
  }
  jackknife <- lapply(unique(subset$pseudotime_bin), function(bin) {
    leave_one_bin_out <- subset[pseudotime_bin != bin]
    vapply(names(figure5_axis_score_columns), function(axis_name) {
      patient_axis_time(leave_one_bin_out, figure5_axis_score_columns[[axis_name]], axis_name)
    }, numeric(1))
  })
  jackknife <- if (length(jackknife)) do.call(rbind, jackknife) else matrix(numeric(), nrow = 0L, ncol = length(figure5_axis_score_columns))
  for (comparison in names(comparison_map)) {
    pair <- comparison_map[[comparison]]
    delta <- times[[pair[[2L]]]] - times[[pair[[1L]]]]
    jackknife_delta <- if (nrow(jackknife)) jackknife[, pair[[2L]]] - jackknife[, pair[[1L]]] else numeric()
    jackknife_delta <- jackknife_delta[is.finite(jackknife_delta)]
    n_jackknife <- length(jackknife_delta)
    se <- if (n_jackknife >= 3L) stats::sd(jackknife_delta) * (n_jackknife - 1) / sqrt(n_jackknife) else NA_real_
    if (is.finite(se)) se <- max(se, 0.005)
    patient_rows[[paste(patient, comparison)]] <- data.table(
      patient_id = patient,
      patient_id_source = paste(unique(subset$patient_id_source), collapse = ";"),
      analysis_level = "evaluable patient-token-level",
      dataset_id = unique(subset$dataset_id)[[1L]],
      comparison = comparison,
      upstream_time = times[[pair[[1L]]]],
      downstream_time = times[[pair[[2L]]]],
      delta = delta,
      se = se,
      ci_lower = delta - 1.96 * se,
      ci_upper = delta + 1.96 * se,
      n_cells = sum(subset$n_cells),
      n_bins = n_bins,
      n_jackknife_valid = n_jackknife,
      coverage_min = coverage[[1L]],
      coverage_max = coverage[[2L]],
      coverage_ok = coverage_ok,
      same_direction = is.finite(delta) && delta > figure5_temporal_parameters$precedence_tolerance,
      se_method = "leave-one-pseudotime-bin-out jackknife: sd(theta_leave_one_out)*(m-1)/sqrt(m), with 0.005 grid floor"
    )
  }
}
patients <- rbindlist(patient_rows, fill = TRUE)

run_meta <- function(subset, stratum = "overall") {
  usable <- subset[is.finite(delta) & is.finite(se) & se > 0]
  comparison <- unique(subset$comparison)[[1L]]
  if (nrow(usable) < 2L) {
    return(data.table(comparison = comparison, stratum = stratum, status = "Not available",
                      n_evaluable_tokens = nrow(usable)))
  }
  fit <- metafor::rma(yi = delta, sei = se, data = usable, method = "REML")
  data.table(
    comparison = comparison,
    stratum = stratum,
    status = "fitted",
    n_evaluable_tokens = nrow(usable),
    pooled_delta = as.numeric(fit$b),
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
    p_value = fit$pval,
    I2 = fit$I2,
    tau2 = fit$tau2,
    same_direction_fraction = mean(usable$delta > figure5_temporal_parameters$precedence_tolerance)
  )
}

meta_rows <- list()
for (comparison in names(comparison_map)) {
  comparison_name <- comparison
  subset <- patients[comparison == comparison_name]
  meta_rows[[comparison]] <- run_meta(subset)
  for (dataset in unique(subset$dataset_id)) {
    if (nrow(subset[dataset_id == dataset]) >= 2L) {
      meta_rows[[paste(comparison, dataset)]] <- run_meta(subset[dataset_id == dataset], paste0("dataset:", dataset))
    }
  }
}
meta <- rbindlist(meta_rows, fill = TRUE)
figure5_write_tsv(patients, file.path(paths$metadata, "figure5g_patient_temporal_differences.tsv"))
figure5_write_tsv(meta, file.path(paths$metadata, "figure5g_patient_meta_analysis.tsv"))
figure5_write_json(list(
  analysis_level = "evaluable patient-token-level; identifiers are derived from sample-name patient tokens and are not an independent clinical validation cohort",
  patient_id_sources = unique(patients$patient_id_source),
  model = "metafor::rma random-effects REML",
  axis_model = "unified temporal-landmark GAM with REML and corrected maximum-slope definition",
  minimum_pseudotime_bins = 5L,
  coverage_rule = "minimum <=0.10 and maximum >=0.90",
  delta_definition = "downstream corrected maximum-slope time minus upstream corrected maximum-slope time; delta >0.005 indicates upstream earlier",
  meta_results = split(meta, seq_len(nrow(meta))),
  packages = figure5_package_versions(c("metafor", "mgcv"))
), file.path(paths$metadata, "figure5g_patient_temporal_report.json"))
