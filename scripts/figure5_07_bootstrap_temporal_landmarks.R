#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
cli <- figure5_cli(list(n_bootstrap = 1000L, seed = 20260805L, n_workers = 4L))
if (!is.finite(cli$n_bootstrap) || cli$n_bootstrap < 1L) stop("n_bootstrap must be positive")
if (!is.finite(cli$seed)) stop("seed must be finite")

pseudo <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
selection <- select_figure5_primary_pseudobulk(pseudo)
primary <- selection$data
unit <- selection$unit
if (uniqueN(primary$patient_id) < 2L) stop("Fewer than two eligible patient/sample units for bootstrap")

resample_stratified <- function(dt) {
  rbindlist(lapply(unique(dt$dataset_id), function(ds) {
    subset <- dt[dataset_id == ds]
    ids <- unique(subset$patient_id)
    sampled <- sample(ids, length(ids), replace = TRUE)
    rbindlist(lapply(seq_along(sampled), function(j) {
      copy_row <- copy(subset[patient_id == sampled[[j]]])
      copy_row[, patient_id := paste0(sampled[[j]], "::bootstrap_copy_", j)]
      copy_row
    }), fill = TRUE)
  }), fill = TRUE)
}

bootstrap_axis <- function(dt, axis_name) {
  score_col <- figure5_axis_score_columns[[axis_name]]
  fit <- tryCatch(
    landmarks_from_table(dt, score_col, axis_name, k = 5L, adjusted = FALSE, use_random_effect = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    diagnostics <- empty_temporal_diagnostics()
    diagnostics$failure_reason <- "model_exception"
    landmarks <- empty_temporal_landmarks(diagnostics)
    model_fitted <- FALSE
  } else {
    landmarks <- fit$landmarks
    diagnostics <- fit$diagnostics
    model_fitted <- !is.null(fit$model)
  }
  list(landmarks = landmarks, diagnostics = diagnostics, model_fitted = model_fitted)
}

run_iteration <- function(iteration) {
  iteration_seed <- as.integer(cli$seed + iteration)
  set.seed(iteration_seed)
  bootstrap_data <- resample_stratified(primary)
  landmark_rows <- vector("list", length(figure5_axis_score_columns))
  diagnostic_rows <- vector("list", length(figure5_axis_score_columns))
  for (position in seq_along(figure5_axis_score_columns)) {
    axis_name <- names(figure5_axis_score_columns)[[position]]
    result <- bootstrap_axis(bootstrap_data, axis_name)
    landmarks <- result$landmarks
    diagnostics <- result$diagnostics
    landmark_rows[[position]] <- data.table(
      iteration = iteration,
      iteration_seed = iteration_seed,
      bootstrap_unit = unit,
      axis = axis_name,
      landmark = names(landmarks),
      time = as.numeric(unlist(landmarks)),
      n_units = uniqueN(bootstrap_data$patient_id),
      n_datasets = uniqueN(bootstrap_data$dataset_id),
      base_seed = cli$seed,
      gam_method = "REML"
    )
    diagnostic_rows[[position]] <- as.data.table(diagnostics)[, `:=`(
      iteration = iteration,
      iteration_seed = iteration_seed,
      bootstrap_unit = unit,
      axis = axis_name,
      n_units = uniqueN(bootstrap_data$patient_id),
      n_datasets = uniqueN(bootstrap_data$dataset_id),
      n_rows = nrow(bootstrap_data),
      n_unique_pseudotime = uniqueN(bootstrap_data$pseudotime),
      model_fitted = result$model_fitted,
      boundary_hit = isTRUE(diagnostics$maximum_slope_boundary_hit) ||
        (is.finite(result$landmarks$maximum_slope_time) &&
           result$landmarks$maximum_slope_time <= figure5_temporal_parameters$baseline_end)
    )]
  }
  list(landmarks = rbindlist(landmark_rows, fill = TRUE),
       diagnostics = rbindlist(diagnostic_rows, fill = TRUE))
}

set.seed(cli$seed)
n_workers <- min(max(1L, as.integer(cli$n_workers)), max(1L, parallel::detectCores(logical = FALSE) - 1L))
if (n_workers > 1L && cli$n_bootstrap > 1L) {
  cluster <- parallel::makeCluster(n_workers)
  parallel::clusterExport(cluster, c("root", "cli", "primary", "unit", "resample_stratified",
                                     "bootstrap_axis", "run_iteration"), envir = environment())
  parallel::clusterEvalQ(cluster, {
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(mgcv))
    source(file.path(root, "scripts", "figure5_temporal_core.R"))
    NULL
  })
  iteration_results <- tryCatch(
    parallel::parLapplyLB(cluster, seq_len(cli$n_bootstrap), run_iteration),
    finally = parallel::stopCluster(cluster)
  )
} else {
  iteration_results <- lapply(seq_len(cli$n_bootstrap), run_iteration)
}

boot <- rbindlist(lapply(iteration_results, `[[`, "landmarks"), fill = TRUE)
diagnostics <- rbindlist(lapply(iteration_results, `[[`, "diagnostics"), fill = TRUE)
finite_summary <- boot[, .(
  n_total = .N,
  n_finite = sum(is.finite(time)),
  finite_fraction = mean(is.finite(time))
), by = .(axis, landmark)]
onset_summary <- boot[landmark == "onset_time", {
  finite <- time[is.finite(time)]
  list(
    finite_fraction = mean(is.finite(time)),
    onset_na_fraction = mean(!is.finite(time)),
    median_onset = if (length(finite)) stats::median(finite) else NA_real_,
    onset_q025 = if (length(finite)) unname(stats::quantile(finite, 0.025, type = 8)) else NA_real_,
    onset_q975 = if (length(finite)) unname(stats::quantile(finite, 0.975, type = 8)) else NA_real_
  )
}, by = axis]
onset_summary[, onset_stability := ifelse(finite_fraction >= 0.80, "Stable", "Unstable")]
diagnostic_summary <- diagnostics[, .(
  coverage_failure_fraction = mean(!coverage_ok),
  boundary_hit_fraction = mean(boundary_hit),
  model_failure_fraction = mean(!model_fitted),
  onset_found_fraction = mean(onset_found)
), by = axis]

development_suffix <- if (cli$n_bootstrap >= 1000L) "" else sprintf("_development_%d", cli$n_bootstrap)
landmark_name <- paste0("figure5_bootstrap_temporal_landmarks", development_suffix)
diagnostic_name <- paste0("figure5_bootstrap_temporal_landmark_diagnostics", development_suffix)
figure5_write_tsv(boot, file.path(paths$metadata, paste0(landmark_name, ".tsv.gz")))
saveRDS(boot, file.path(paths$processed, paste0(landmark_name, ".rds")))
figure5_write_tsv(diagnostics, file.path(paths$metadata, paste0(diagnostic_name, ".tsv.gz")))
figure5_write_json(list(
  n_bootstrap_requested = cli$n_bootstrap,
  n_bootstrap_complete = uniqueN(boot$iteration),
  bootstrap_unit = unit,
  analysis_unit = selection$analysis_unit_note,
  coverage_fallback_reason = selection$fallback_reason,
  stratified_by = "dataset_id",
  bootstrap_model = "same fit_landmark_gam REML structure as main analysis: cubic regression spline, dataset fixed effect, estimable patient/sample random effect",
  resampling_cells_as_independent_units = FALSE,
  coverage_rule = "observed minimum <=0.10 and observed maximum >=0.90; failures return NA without extrapolation",
  random_seed = cli$seed,
  iteration_seed_rule = "base_seed + iteration",
  n_workers = n_workers,
  fixed_parameters = figure5_temporal_parameters,
  finite_fraction = split(finite_summary, seq_len(nrow(finite_summary))),
  onset_summary = split(onset_summary, seq_len(nrow(onset_summary))),
  diagnostic_summary = split(diagnostic_summary, seq_len(nrow(diagnostic_summary))),
  review_risk = if (cli$n_bootstrap < 1000L) "Development run below 1,000 replicates" else
    onset_summary[onset_stability == "Unstable", paste(axis, "onset finite fraction <0.80; onset excluded from main-figure claims")]
), file.path(paths$metadata, paste0(landmark_name, "_report.json")))

message(sprintf("Figure 5 corrected bootstrap complete: %d iterations, %d workers, %s unit.",
                cli$n_bootstrap, n_workers, unit))
