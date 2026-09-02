#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
pseudo <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
boot <- as.data.table(readRDS(file.path(paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))

pair_probabilities <- function(times, method_name, evidence_unit,
                               independence_status = "independent sensitivity view",
                               reported_n_units = NULL, unit_type = "fitted units") {
  if (is.null(reported_n_units)) reported_n_units <- nrow(times)
  out <- lapply(seq_len(nrow(precedence_pairs)), function(i) {
    pair <- precedence_pairs[i]
    statistics <- tie_aware_precedence_probability(times[[pair$upstream_axis]], times[[pair$downstream_axis]])
    grade <- classify_precedence_evidence(
      statistics$probability, statistics$valid_fraction, statistics$tie_fraction,
      statistics$delta_q025, statistics$delta_q975
    )
    weight <- if (grepl("non-independent", independence_status, fixed = TRUE)) {
      "non-independent sensitivity view"
    } else if (statistics$n_valid < 10L) {
      "descriptive only; small n_units"
    } else {
      "supporting sensitivity analysis"
    }
    data.table(
      method = method_name,
      comparison = pair$comparison,
      probability = statistics$probability,
      n_units = as.integer(reported_n_units),
      n_valid = statistics$n_valid,
      n_tied = statistics$n_tied,
      tie_fraction = statistics$tie_fraction,
      median_delta = statistics$median_delta,
      delta_q025 = statistics$delta_q025,
      delta_q975 = statistics$delta_q975,
      status = grade,
      evidence_unit = evidence_unit,
      unit_type = unit_type,
      independence_status = independence_status,
      evidence_weight = weight
    )
  })
  rbindlist(out)
}

unit_landmarks <- function(dt, unit_col = "patient_id", score_columns = figure5_axis_score_columns, k = 5L) {
  unit_axis_time <- function(data, score_col, axis_name, requested_k) {
    work <- data[is.finite(pseudotime) & is.finite(get(score_col))]
    coverage <- if (nrow(work)) range(work$pseudotime, na.rm = TRUE) else c(NA_real_, NA_real_)
    if (nrow(work) < 5L || uniqueN(work$pseudotime) < 5L ||
        !all(is.finite(coverage)) || coverage[[1L]] > figure5_temporal_parameters$baseline_end || coverage[[2L]] < 0.90) {
      return(NA_real_)
    }
    fit <- landmarks_from_table(work, score_col, axis_name, k = requested_k,
                                adjusted = FALSE, use_random_effect = FALSE)
    fit$landmarks$maximum_slope_time
  }
  units <- unique(dt[[unit_col]])
  rows <- lapply(units, function(unit_id) {
    subset <- dt[get(unit_col) == unit_id]
    if (nrow(subset) < 5L || uniqueN(subset$pseudotime) < 5L) return(NULL)
    times <- vapply(names(score_columns), function(axis_name) {
      score_col <- score_columns[[axis_name]]
      if (!score_col %in% names(subset) || sum(is.finite(subset[[score_col]])) < 5L) return(NA_real_)
      unit_axis_time(subset, score_col, axis_name, k)
    }, numeric(1))
    as.data.table(as.list(times))[, unit := as.character(unit_id)]
  })
  rbindlist(rows, fill = TRUE)
}

selection <- select_figure5_primary_pseudobulk(pseudo)
primary <- selection$data
sample_dt <- pseudo[aggregation_unit == "sample" & method == "main/consensus pseudotime" & patient_meta_eligible == TRUE]
concordance <- list()
main_wide <- dcast(boot[landmark == "maximum_slope_time"], iteration ~ axis, value.var = "time")
concordance[["main"]] <- pair_probabilities(
  main_wide, "main/consensus pseudotime", "corrected stratified bootstrap replicates",
  independence_status = "primary corrected bootstrap",
  reported_n_units = uniqueN(primary$patient_id), unit_type = paste0(selection$unit, " tokens resampled within dataset")
)
concordance[["patient"]] <- pair_probabilities(
  main_wide, "patient-pseudobulk", "same corrected stratified bootstrap replicates",
  independence_status = "non-independent duplicate of main/consensus corrected bootstrap",
  reported_n_units = uniqueN(primary$patient_id), unit_type = paste0(selection$unit, " tokens resampled within dataset")
)

for (method_name in c("Monocle3", "Slingshot scanVI", "Slingshot hepatocyte PCA", "CytoTRACE2")) {
  subset <- pseudo[aggregation_unit == "patient" & method == method_name & patient_meta_eligible == TRUE]
  times <- unit_landmarks(subset)
  concordance[[method_name]] <- pair_probabilities(
    times, method_name, "evaluable patient-token-level REML fits",
    reported_n_units = nrow(times), unit_type = "evaluable patient tokens"
  )
}

sample_times <- unit_landmarks(sample_dt, unit_col = "sample_id")
concordance[["sample"]] <- pair_probabilities(
  sample_times, "sample-balanced", "evaluable sample-token-level REML fits",
  reported_n_units = nrow(sample_times), unit_type = "evaluable sample tokens"
)

scenario_specs <- list(
  "CNV-strict" = list(data = primary[cnv_strict_fraction >= 0.5 | pseudotime <= 0.2], scores = figure5_axis_score_columns),
  "no-proliferation" = list(data = primary[!is.finite(proliferation_score) | proliferation_score <= quantile(proliferation_score, 0.90, na.rm = TRUE)], scores = figure5_axis_score_columns),
  "no-generic-stress" = list(data = primary, scores = c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score_no_generic", sox4_stabilization = "sox4_stabilization_score")),
  "GAM k=4" = list(data = primary, scores = figure5_axis_score_columns, k = 4L),
  "GAM k=6" = list(data = primary, scores = figure5_axis_score_columns, k = 6L)
)
for (name in names(scenario_specs)) {
  specification <- scenario_specs[[name]]
  times <- unit_landmarks(specification$data, score_columns = specification$scores,
                          k = if (is.null(specification$k)) 5L else specification$k)
  concordance[[name]] <- pair_probabilities(
    times, name, "evaluable patient-token-level REML fits",
    reported_n_units = nrow(times), unit_type = "evaluable patient tokens"
  )
}

omit_prob <- function(dt, omit_col, label) {
  rows <- lapply(unique(dt[[omit_col]]), function(omit) {
    subset <- dt[get(omit_col) != omit]
    times <- vapply(names(figure5_axis_score_columns), function(axis_name) {
      landmarks_from_table(subset, figure5_axis_score_columns[[axis_name]], axis_name,
                           use_random_effect = TRUE)$landmarks$maximum_slope_time
    }, numeric(1))
    as.data.table(as.list(times))[, omitted := as.character(omit)]
  })
  times <- rbindlist(rows, fill = TRUE)
  pair_probabilities(times, label, paste0("leave-one-", omit_col, "-out REML fits"),
                     reported_n_units = nrow(times), unit_type = paste0("omitted ", omit_col, " scenarios"))
}
concordance[["LODO"]] <- omit_prob(primary, "dataset_id", "leave-one-dataset-out")
concordance[["LOSO"]] <- omit_prob(sample_dt, "sample_id", "leave-one-sample-out")

unavailable <- rbindlist(lapply(c("DPT", "CellRank pseudotime"), function(method_name) data.table(
  method = method_name,
  comparison = precedence_pairs$comparison,
  probability = NA_real_, n_units = 0L, n_valid = 0L, n_tied = 0L,
  tie_fraction = NA_real_, median_delta = NA_real_, delta_q025 = NA_real_, delta_q975 = NA_real_,
  status = "Not available", evidence_unit = "method unavailable", unit_type = "unavailable",
  independence_status = "not evaluated", evidence_weight = "none"
)))

result <- rbindlist(c(concordance, list(unavailable)), fill = TRUE)
result[, landmark := "maximum_slope"]
figure5_write_tsv(result, file.path(paths$metadata, "figure5f_method_concordance.tsv"))
figure5_write_tsv(result, file.path(paths$metadata, "figure5_sensitivity_summary.tsv"))
figure5_write_json(list(
  methods_tested = unique(result$method),
  unavailable_methods = c("DPT", "CellRank pseudotime"),
  gam_method = "REML",
  minimum_distinct_pseudotime_positions = 5L,
  coverage_rule = "minimum <=0.10 and maximum >=0.90",
  precedence_rule = "tie-aware with tolerance 0.005",
  non_independence = "main/consensus and patient-pseudobulk rows are two labels for the same corrected bootstrap result",
  small_n_rule = "fewer than 10 valid fitted units is descriptive only",
  primary_landmark = "maximum_slope"
), file.path(paths$metadata, "figure5f_method_concordance_report.json"))
