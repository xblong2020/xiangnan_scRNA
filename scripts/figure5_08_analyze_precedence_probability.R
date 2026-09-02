#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_onset_fix_paths(root)
boot <- as.data.table(readRDS(file.path(paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
if (uniqueN(boot$iteration) < 1000L) stop("Formal precedence analysis requires the corrected 1,000-iteration bootstrap")

landmark_map <- c(onset = "onset_time", maximum_slope = "maximum_slope_time", t50 = "t50")
onset_finite <- boot[landmark == "onset_time", .(onset_finite_fraction = mean(is.finite(time))), by = axis]
onset_fraction <- setNames(onset_finite$onset_finite_fraction, onset_finite$axis)
n_iterations <- uniqueN(boot$iteration)
rows <- list()

for (label in names(landmark_map)) {
  wide <- dcast(boot[landmark == landmark_map[[label]]], iteration ~ axis, value.var = "time")
  for (i in seq_len(nrow(precedence_pairs))) {
    pair <- precedence_pairs[i]
    statistics <- tie_aware_precedence_probability(
      wide[[pair$upstream_axis]], wide[[pair$downstream_axis]],
      tolerance = figure5_temporal_parameters$precedence_tolerance
    )
    onset_pair_fraction <- min(onset_fraction[[pair$upstream_axis]], onset_fraction[[pair$downstream_axis]], na.rm = TRUE)
    onset_for_grade <- if (label == "onset") onset_pair_fraction else NA_real_
    grade <- classify_precedence_evidence(
      statistics$probability, statistics$valid_fraction, statistics$tie_fraction,
      statistics$delta_q025, statistics$delta_q975,
      onset_finite_fraction = onset_for_grade
    )
    mc_lower <- if (is.finite(statistics$Monte_Carlo_SE)) max(0, statistics$probability - 1.96 * statistics$Monte_Carlo_SE) else NA_real_
    mc_upper <- if (is.finite(statistics$Monte_Carlo_SE)) min(1, statistics$probability + 1.96 * statistics$Monte_Carlo_SE) else NA_real_
    rows[[paste(label, pair$comparison)]] <- data.table(
      landmark = label,
      comparison = pair$comparison,
      upstream_axis = pair$upstream_axis,
      downstream_axis = pair$downstream_axis,
      n_bootstrap_total = n_iterations,
      n_valid = statistics$n_valid,
      n_earlier = statistics$n_earlier,
      n_tied = statistics$n_tied,
      n_later = statistics$n_later,
      valid_fraction = statistics$valid_fraction,
      tie_fraction = statistics$tie_fraction,
      probability = statistics$probability,
      median_delta = statistics$median_delta,
      delta_q025 = statistics$delta_q025,
      delta_q975 = statistics$delta_q975,
      Monte_Carlo_SE = statistics$Monte_Carlo_SE,
      mc_interval_lower = mc_lower,
      mc_interval_upper = mc_upper,
      mc_interval_label = "Monte Carlo interval only",
      onset_finite_fraction_pair = onset_pair_fraction,
      evidence_grade = grade,
      is_primary = label == "maximum_slope"
    )
  }
}

precedence <- rbindlist(rows)
figure5_write_tsv(precedence, file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))
primary <- precedence[is_primary == TRUE]
conclusion <- figure5_result_text(primary)
figure5_write_json(list(
  primary_landmark = "maximum_slope",
  tolerance = figure5_temporal_parameters$precedence_tolerance,
  tie_score = 0.5,
  probability_definition = "mean of 1 for delta>tolerance, 0.5 for ties, and 0 for delta<−tolerance",
  delta_definition = "downstream_time - upstream_time",
  intervals = "delta_q025 and delta_q975 are bootstrap distribution quantiles; probability intervals, if shown, are Monte Carlo intervals only",
  evidence_grade_considers = c("valid bootstrap fraction", "tie fraction", "delta distribution", "onset stability for onset comparisons"),
  probabilities = split(precedence, seq_len(nrow(precedence))),
  conclusion = conclusion,
  interpretation_constraint = "Cross-sectional pseudotime supports relative positioning only; it does not establish physical time or causality."
), file.path(paths$metadata, "figure5e_precedence_report.json"))
message(conclusion)
