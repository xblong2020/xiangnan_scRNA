#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

derive_figure5_patient_id <- function(dataset_id, cnv_sample) {
  dataset_id <- as.character(dataset_id)
  cnv_sample <- as.character(cnv_sample)
  patient_token <- cnv_sample
  source <- rep("aggregate_sample_proxy", length(cnv_sample))

  hcc <- dataset_id == "GSE149614" & grepl("^HCC[0-9]+", cnv_sample)
  patient_token[hcc] <- sub("^(HCC[0-9]+).*$", "\\1", cnv_sample[hcc])
  source[hcc] <- "sample_name_patient_token"

  c_patient <- dataset_id == "GSE185477" & grepl("^C[0-9]+", cnv_sample)
  patient_token[c_patient] <- sub("^(C[0-9]+).*$", "\\1", cnv_sample[c_patient])
  source[c_patient] <- "sample_name_patient_token"

  s_patient <- dataset_id == "GSE212046" & grepl("^S[0-9]+$", cnv_sample)
  patient_token[s_patient] <- cnv_sample[s_patient]
  source[s_patient] <- "sample_name_patient_token"

  data.table(
    patient_id = paste(dataset_id, patient_token, sep = "::"),
    patient_id_source = source,
    patient_meta_eligible = source == "sample_name_patient_token"
  )
}

robust_z_by_group <- function(values, groups) {
  x <- suppressWarnings(as.numeric(values))
  groups <- as.character(groups)
  out <- rep(NA_real_, length(x))
  for (group in unique(groups)) {
    idx <- which(groups == group & is.finite(x))
    if (!length(idx)) next
    center <- stats::median(x[idx], na.rm = TRUE)
    spread <- stats::mad(x[idx], center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(spread) || spread <= .Machine$double.eps) {
      out[idx] <- 0
    } else {
      out[idx] <- pmax(-5, pmin(5, (x[idx] - center) / spread))
    }
  }
  out
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

decide_pseudotime_orientation <- function(
  pseudotime,
  mature_identity_score,
  cnv_score,
  malignant_fate,
  malignant_indicator
) {
  pt <- suppressWarnings(as.numeric(pseudotime))
  if (all(!is.finite(pt))) {
    return(list(flipped = FALSE, original_score = NA_real_, oriented_score = NA_real_, evidence = numeric()))
  }
  range_pt <- range(pt, na.rm = TRUE)
  if (diff(range_pt) > 0) pt <- (pt - range_pt[[1]]) / diff(range_pt)
  evidence <- c(
    identity_loss = safe_spearman(pt, -as.numeric(mature_identity_score)),
    cnv = safe_spearman(pt, as.numeric(cnv_score)),
    malignant_fate = safe_spearman(pt, as.numeric(malignant_fate)),
    malignant_fraction = safe_spearman(pt, as.numeric(malignant_indicator))
  )
  score <- mean(evidence[is.finite(evidence)], na.rm = TRUE)
  if (!is.finite(score)) score <- 0
  flipped <- score < 0
  list(
    flipped = flipped,
    original_score = score,
    oriented_score = if (flipped) -score else score,
    pseudotime_oriented = if (flipped) 1 - pt else pt,
    evidence = evidence
  )
}

figure5_temporal_parameters <- list(
  baseline_start = 0,
  baseline_end = 0.10,
  search_start = 0.10,
  slope_search_end = 0.95,
  grid_n = 201L,
  onset_min_run = 5L,
  crossing_min_run = 3L,
  min_effect_fraction = 0.05,
  min_effect_baseline_sd = 0.25,
  derivative_fraction = 0.10,
  precedence_tolerance = 0.005,
  fitted_scale_fallback_fraction = 0.10
)

temporal_landmark_names <- c(
  "onset_time", "t10", "t50", "maximum_slope_time", "extremum_time",
  "peak_time", "plateau_time", "decline_onset"
)

empty_temporal_diagnostics <- function() {
  list(
    baseline_start = figure5_temporal_parameters$baseline_start,
    baseline_end = figure5_temporal_parameters$baseline_end,
    search_start = figure5_temporal_parameters$search_start,
    baseline_value = NA_real_,
    baseline_scale = NA_real_,
    baseline_scale_source = "unavailable",
    peak_time = NA_real_,
    peak_value = NA_real_,
    total_rise = NA_real_,
    effect_threshold = NA_real_,
    maximum_positive_derivative = NA_real_,
    derivative_threshold = NA_real_,
    onset_candidate_count = 0L,
    onset_run_length = 0L,
    onset_found = FALSE,
    maximum_slope_boundary_hit = FALSE,
    curve_coverage_min = NA_real_,
    curve_coverage_max = NA_real_,
    coverage_ok = FALSE,
    failure_reason = "not_evaluated"
  )
}

empty_temporal_landmarks <- function(diagnostics = empty_temporal_diagnostics()) {
  result <- as.list(stats::setNames(rep(NA_real_, length(temporal_landmark_names)), temporal_landmark_names))
  attr(result, "diagnostics") <- diagnostics
  result
}

central_derivative <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  out <- rep(NA_real_, length(x))
  if (length(x) < 3L || length(y) != length(x)) return(out)
  interior <- 2L:(length(x) - 1L)
  denominator <- x[interior + 1L] - x[interior - 1L]
  valid <- is.finite(denominator) & denominator > .Machine$double.eps &
    is.finite(y[interior + 1L]) & is.finite(y[interior - 1L])
  out[interior[valid]] <- (y[interior[valid] + 1L] - y[interior[valid] - 1L]) / denominator[valid]
  out
}

first_persistent_run <- function(condition, min_run = 1L) {
  condition <- as.logical(condition)
  condition[is.na(condition)] <- FALSE
  min_run <- max(1L, as.integer(min_run))
  if (!length(condition)) return(NA_integer_)
  runs <- rle(condition)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  valid <- which(runs$values & runs$lengths >= min_run)
  if (!length(valid)) return(NA_integer_)
  answer <- starts[valid[[1L]]]
  attr(answer, "run_length") <- as.integer(runs$lengths[valid[[1L]]])
  answer
}

persistent_crossing_time <- function(x, y, target, min_run = 3L, direction = 1, eligible = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (is.null(eligible)) eligible <- rep(TRUE, length(x))
  condition <- as.logical(eligible) & is.finite(y) & (sign(direction) * (y - target) >= 0)
  index <- first_persistent_run(condition, min_run = min_run)
  if (is.na(index)) NA_real_ else x[[index]]
}

figure5_robust_scale <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  if (length(x) < 2L) return(NA_real_)
  center <- stats::median(x)
  value <- stats::mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (is.finite(value) && value > .Machine$double.eps) return(value)
  value <- stats::sd(x, na.rm = TRUE)
  if (is.finite(value) && value > .Machine$double.eps) value else NA_real_
}

figure5h_robust_scale <- figure5_robust_scale

resolve_temporal_baseline_scale <- function(raw_values = numeric(), fitted = numeric(),
                                            fallback_fraction = figure5_temporal_parameters$fitted_scale_fallback_fraction) {
  raw_values <- as.numeric(raw_values[is.finite(raw_values)])
  if (length(raw_values) >= 2L) {
    center <- stats::median(raw_values)
    mad_value <- stats::mad(raw_values, center = center, constant = 1.4826, na.rm = TRUE)
    if (is.finite(mad_value) && mad_value > .Machine$double.eps) {
      return(list(scale = mad_value, source = "raw_baseline_mad", n = length(raw_values)))
    }
    sd_value <- stats::sd(raw_values, na.rm = TRUE)
    if (is.finite(sd_value) && sd_value > .Machine$double.eps) {
      return(list(scale = sd_value, source = "raw_baseline_sd", n = length(raw_values)))
    }
  }
  fitted_scale <- figure5_robust_scale(fitted)
  if (is.finite(fitted_scale) && fitted_scale > .Machine$double.eps) {
    return(list(scale = fallback_fraction * fitted_scale,
                source = sprintf("fitted_curve_robust_scale_x%.2f", fallback_fraction), n = length(raw_values)))
  }
  list(scale = NA_real_, source = "unavailable", n = length(raw_values))
}

estimate_temporal_landmarks <- function(
  pseudotime,
  fitted,
  axis = NULL,
  baseline_scale = NA_real_,
  baseline_scale_source = if (is.finite(baseline_scale)) "provided" else "unavailable",
  observed_range = NULL,
  parameters = figure5_temporal_parameters
) {
  keep <- is.finite(pseudotime) & is.finite(fitted)
  curve <- data.table(pseudotime = as.numeric(pseudotime[keep]), fitted = as.numeric(fitted[keep]))
  if (nrow(curve)) curve <- curve[, .(fitted = stats::median(fitted)), by = pseudotime][order(pseudotime)]
  x <- curve$pseudotime
  y <- curve$fitted
  diagnostics <- empty_temporal_diagnostics()
  diagnostics$baseline_start <- parameters$baseline_start
  diagnostics$baseline_end <- parameters$baseline_end
  diagnostics$search_start <- parameters$search_start
  diagnostics$baseline_scale <- if (is.finite(baseline_scale)) as.numeric(baseline_scale) else NA_real_
  diagnostics$baseline_scale_source <- baseline_scale_source

  coverage <- if (is.null(observed_range)) range(x, na.rm = TRUE) else as.numeric(observed_range)
  if (length(coverage) >= 2L && all(is.finite(coverage[1:2]))) {
    diagnostics$curve_coverage_min <- min(coverage[1:2])
    diagnostics$curve_coverage_max <- max(coverage[1:2])
  }
  diagnostics$coverage_ok <- is.finite(diagnostics$curve_coverage_min) &&
    is.finite(diagnostics$curve_coverage_max) &&
    diagnostics$curve_coverage_min <= parameters$baseline_end + 1e-8 &&
    diagnostics$curve_coverage_max >= 0.90 - 1e-8
  if (!isTRUE(diagnostics$coverage_ok)) {
    early <- !is.finite(diagnostics$curve_coverage_min) || diagnostics$curve_coverage_min > parameters$baseline_end + 1e-8
    late <- !is.finite(diagnostics$curve_coverage_max) || diagnostics$curve_coverage_max < 0.90 - 1e-8
    diagnostics$failure_reason <- if (early && late) "insufficient_early_and_late_coverage" else if (early) "insufficient_early_coverage" else "insufficient_late_coverage"
    return(empty_temporal_landmarks(diagnostics))
  }
  if (length(x) < max(5L, parameters$onset_min_run) || diff(range(x)) <= .Machine$double.eps) {
    diagnostics$failure_reason <- "insufficient_curve_grid"
    return(empty_temporal_landmarks(diagnostics))
  }

  baseline_index <- which(x >= parameters$baseline_start - 1e-8 & x <= parameters$baseline_end + 1e-8)
  if (length(baseline_index) < 2L) {
    diagnostics$failure_reason <- "insufficient_baseline_grid"
    return(empty_temporal_landmarks(diagnostics))
  }
  diagnostics$baseline_value <- stats::median(y[baseline_index], na.rm = TRUE)
  if (!is.finite(diagnostics$baseline_scale)) {
    fallback <- resolve_temporal_baseline_scale(y[baseline_index], y, parameters$fitted_scale_fallback_fraction)
    diagnostics$baseline_scale <- fallback$scale
    diagnostics$baseline_scale_source <- paste0("fitted_baseline_", fallback$source)
  }

  post_baseline <- which(x > parameters$search_start + 1e-8)
  if (!length(post_baseline)) {
    diagnostics$failure_reason <- "no_post_baseline_grid"
    return(empty_temporal_landmarks(diagnostics))
  }
  peak_index <- post_baseline[[which.max(y[post_baseline])]]
  diagnostics$peak_time <- x[[peak_index]]
  diagnostics$peak_value <- y[[peak_index]]
  diagnostics$total_rise <- diagnostics$peak_value - diagnostics$baseline_value
  if (!is.finite(diagnostics$total_rise) || diagnostics$total_rise <= .Machine$double.eps) {
    diagnostics$failure_reason <- "no_positive_baseline_to_peak_rise"
    return(empty_temporal_landmarks(diagnostics))
  }
  if (!is.finite(diagnostics$baseline_scale) || diagnostics$baseline_scale <= 0) {
    diagnostics$failure_reason <- "baseline_scale_unavailable"
    return(empty_temporal_landmarks(diagnostics))
  }

  derivative <- central_derivative(x, y)
  slope_limit <- min(diagnostics$peak_time, parameters$slope_search_end)
  slope_eligible <- x > parameters$search_start + 1e-8 & x <= slope_limit + 1e-8 & is.finite(derivative)
  positive_derivatives <- derivative[slope_eligible & derivative > 0]
  if (length(positive_derivatives)) diagnostics$maximum_positive_derivative <- max(positive_derivatives)
  if (is.finite(diagnostics$maximum_positive_derivative)) {
    diagnostics$derivative_threshold <- parameters$derivative_fraction * diagnostics$maximum_positive_derivative
  }
  diagnostics$effect_threshold <- max(
    parameters$min_effect_fraction * diagnostics$total_rise,
    parameters$min_effect_baseline_sd * diagnostics$baseline_scale
  )

  onset_condition <- slope_eligible &
    (y - diagnostics$baseline_value >= diagnostics$effect_threshold) &
    derivative > 0 & derivative >= diagnostics$derivative_threshold
  onset_condition[is.na(onset_condition)] <- FALSE
  diagnostics$onset_candidate_count <- sum(onset_condition)
  onset_index <- first_persistent_run(onset_condition, parameters$onset_min_run)
  if (!is.na(onset_index)) {
    diagnostics$onset_run_length <- attr(onset_index, "run_length")
    diagnostics$onset_found <- TRUE
  }

  crossing_eligible <- x > parameters$search_start + 1e-8 & x <= diagnostics$peak_time + 1e-8
  t10 <- persistent_crossing_time(
    x, y, diagnostics$baseline_value + 0.10 * diagnostics$total_rise,
    min_run = parameters$crossing_min_run, direction = 1, eligible = crossing_eligible
  )
  t50 <- persistent_crossing_time(
    x, y, diagnostics$baseline_value + 0.50 * diagnostics$total_rise,
    min_run = parameters$crossing_min_run, direction = 1, eligible = crossing_eligible
  )

  maximum_slope_time <- NA_real_
  slope_indices <- which(slope_eligible & derivative > 0)
  if (length(slope_indices)) {
    maximum_slope_index <- slope_indices[[which.max(derivative[slope_indices])]]
    rightmost_slope_index <- max(slope_indices)
    diagnostics$maximum_slope_boundary_hit <- maximum_slope_index == rightmost_slope_index &&
      x[[maximum_slope_index]] >= slope_limit - 1e-8
    if (!diagnostics$maximum_slope_boundary_hit) maximum_slope_time <- x[[maximum_slope_index]]
  }

  plateau_time <- NA_real_
  if (is.finite(maximum_slope_time) && is.finite(diagnostics$derivative_threshold)) {
    plateau_condition <- x > maximum_slope_time + 1e-8 & is.finite(derivative) &
      abs(derivative) <= diagnostics$derivative_threshold
    plateau_index <- first_persistent_run(plateau_condition, parameters$onset_min_run)
    if (!is.na(plateau_index)) plateau_time <- x[[plateau_index]]
  }

  decline_onset <- NA_real_
  if (is.finite(diagnostics$derivative_threshold)) {
    decline_condition <- x > diagnostics$peak_time + 1e-8 & is.finite(derivative) &
      derivative <= -diagnostics$derivative_threshold
    decline_index <- first_persistent_run(decline_condition, parameters$onset_min_run)
    if (!is.na(decline_index)) decline_onset <- x[[decline_index]]
  }

  failure_reasons <- character()
  if (!diagnostics$onset_found) failure_reasons <- c(failure_reasons, "no_sustained_onset")
  if (diagnostics$maximum_slope_boundary_hit) {
    failure_reasons <- c(failure_reasons, "maximum_slope_at_search_boundary")
  } else if (!is.finite(maximum_slope_time)) {
    failure_reasons <- c(failure_reasons, "no_post_baseline_positive_derivative")
  }
  diagnostics$failure_reason <- if (length(failure_reasons)) paste(failure_reasons, collapse = ";") else ""
  result <- list(
    onset_time = if (diagnostics$onset_found) x[[as.integer(onset_index)]] else NA_real_,
    t10 = t10,
    t50 = t50,
    maximum_slope_time = maximum_slope_time,
    extremum_time = diagnostics$peak_time,
    peak_time = diagnostics$peak_time,
    plateau_time = plateau_time,
    decline_onset = decline_onset
  )
  result <- result[temporal_landmark_names]
  attr(result, "diagnostics") <- diagnostics
  result
}

tie_aware_precedence_probability <- function(upstream_time, downstream_time,
                                             tolerance = figure5_temporal_parameters$precedence_tolerance) {
  upstream_time <- as.numeric(upstream_time)
  downstream_time <- as.numeric(downstream_time)
  valid <- is.finite(upstream_time) & is.finite(downstream_time)
  delta <- downstream_time[valid] - upstream_time[valid]
  scores <- ifelse(delta > tolerance, 1, ifelse(abs(delta) <= tolerance, 0.5, 0))
  n_valid <- length(delta)
  list(
    n_valid = as.integer(n_valid),
    n_earlier = as.integer(sum(delta > tolerance)),
    n_tied = as.integer(sum(abs(delta) <= tolerance)),
    n_later = as.integer(sum(delta < -tolerance)),
    tie_fraction = if (n_valid) mean(abs(delta) <= tolerance) else NA_real_,
    probability = if (n_valid) mean(scores) else NA_real_,
    median_delta = if (n_valid) stats::median(delta) else NA_real_,
    delta_q025 = if (n_valid) unname(stats::quantile(delta, 0.025, type = 8)) else NA_real_,
    delta_q975 = if (n_valid) unname(stats::quantile(delta, 0.975, type = 8)) else NA_real_,
    Monte_Carlo_SE = if (n_valid > 1L) stats::sd(scores) / sqrt(n_valid) else NA_real_,
    valid_fraction = if (length(upstream_time)) n_valid / length(upstream_time) else NA_real_
  )
}

# Find the first sustained, directionally supported activation point. This is
# deliberately separate from onset_time: early movement in the opposite
# direction cannot satisfy the high-score activation contract.
find_directional_activation <- function(pseudotime, fitted, derivative,
                                         derivative_lower, derivative_upper,
                                         direction = 1, min_change_sd = 0.25,
                                         min_run = 5L, baseline_fraction = 0.10) {
  keep <- is.finite(pseudotime) & is.finite(fitted) & is.finite(derivative) &
    is.finite(derivative_lower) & is.finite(derivative_upper)
  if (sum(keep) < max(5L, min_run)) {
    return(list(found = FALSE, time = NA_real_, index = NA_integer_,
                baseline = NA_real_, baseline_sd = NA_real_, criteria = logical()))
  }
  x <- as.numeric(pseudotime[keep])
  y <- as.numeric(fitted[keep])
  d <- as.numeric(derivative[keep])
  dl <- as.numeric(derivative_lower[keep])
  du <- as.numeric(derivative_upper[keep])
  ord <- order(x)
  x <- x[ord]; y <- y[ord]; d <- d[ord]; dl <- dl[ord]; du <- du[ord]

  baseline_n <- max(5L, min(length(y), ceiling(length(y) * baseline_fraction)))
  baseline <- stats::median(y[seq_len(baseline_n)], na.rm = TRUE)
  baseline_sd <- figure5h_robust_scale(y)
  direction <- ifelse(direction >= 0, 1, -1)
  derivative_direction_ok <- direction * d > 0
  derivative_ci_ok <- if (direction > 0) dl > 0 else du < 0
  change_ok <- direction * (y - baseline) >= min_change_sd * baseline_sd
  criteria <- derivative_direction_ok & derivative_ci_ok & change_ok
  runs <- rle(criteria)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  valid_run <- which(runs$values & runs$lengths >= min_run)
  if (!length(valid_run)) {
    return(list(found = FALSE, time = NA_real_, index = NA_integer_,
                baseline = baseline, baseline_sd = baseline_sd, criteria = criteria))
  }
  first_index <- starts[valid_run[[1]]]
  list(found = TRUE, time = x[[first_index]], index = first_index,
       baseline = baseline, baseline_sd = baseline_sd, criteria = criteria)
}

resolve_figure5h_start <- function(t10, t10_lower, t10_upper,
                                   t10_finite_fraction, directional = NULL,
                                   onset = NA_real_, onset_lower = NA_real_, onset_upper = NA_real_,
                                   onset_finite_fraction = 0,
                                   max_t10_ci_width = 0.25,
                                   max_onset_ci_width = 0.25,
                                   min_t10_finite_fraction = 0.80,
                                   min_onset_finite_fraction = 0.80) {
  onset_stable <- is.finite(onset) && is.finite(onset_lower) && is.finite(onset_upper) &&
    is.finite(onset_finite_fraction) && onset_finite_fraction >= min_onset_finite_fraction &&
    (onset_upper - onset_lower) <= max_onset_ci_width
  if (onset_stable) {
    return(list(start = onset, start_lower = onset_lower, start_upper = onset_upper,
                method = "corrected_bootstrap_onset", boundary_status = "resolved"))
  }
  t10_stable <- is.finite(t10) && is.finite(t10_lower) && is.finite(t10_upper) &&
    is.finite(t10_finite_fraction) && t10_finite_fraction >= min_t10_finite_fraction &&
    (t10_upper - t10_lower) <= max_t10_ci_width
  directional_found <- !is.null(directional) && isTRUE(directional$found) && is.finite(directional$time)
  if (t10_stable && directional_found && t10 >= directional$time - 1e-8) {
    return(list(start = t10, start_lower = t10_lower, start_upper = t10_upper,
                method = "bootstrap_t10", boundary_status = "resolved"))
  }
  if (directional_found) {
    # When t10 precedes the directionally supported rise, its lower CI is a
    # boundary-sensitive onset estimate rather than a valid activation range.
    # Keep a local uncertainty interval around the directional point instead
    # of extending the formal band back to pseudotime zero.
    fallback_half_width <- if (t10_stable) max(0.05, (t10_upper - t10_lower)) else 0.05
    lower <- max(0, directional$time - fallback_half_width)
    upper <- if (is.finite(t10_upper)) max(t10_upper, directional$time) else directional$time
    return(list(start = directional$time, start_lower = lower, start_upper = upper,
                method = "directional_activation_fallback",
                boundary_status = if (t10_stable) "t10_pre_directional" else "resolved"))
  }
  if (t10_stable) {
    return(list(start = t10, start_lower = t10_lower, start_upper = t10_upper,
                method = "bootstrap_t10", boundary_status = "resolved"))
  }
  list(start = NA_real_, start_lower = NA_real_, start_upper = NA_real_,
       method = "boundary unresolved", boundary_status = "boundary unresolved")
}

# Convert a smoothed score curve into a continuous prominence band. The curve
# fades at the right edge instead of assigning a discrete programme end.
build_figure5h_activity_profile <- function(pseudotime, fitted, start, t50,
                                             maximum_slope, decline_onset = NA_real_, direction = 1,
                                             baseline_fraction = 0.10) {
  keep <- is.finite(pseudotime) & is.finite(fitted)
  x <- as.numeric(pseudotime[keep]); y <- as.numeric(fitted[keep])
  ord <- order(x); x <- x[ord]; y <- y[ord]
  if (!length(x)) return(data.table())
  baseline_n <- max(5L, min(length(y), ceiling(length(y) * baseline_fraction)))
  baseline <- stats::median(y[seq_len(baseline_n)], na.rm = TRUE)
  scale_value <- figure5h_robust_scale(y)
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  direction <- ifelse(direction >= 0, 1, -1)
  positive_change <- pmax(0, direction * (y - baseline))
  high_reference <- stats::quantile(positive_change[positive_change > 0], 0.90, na.rm = TRUE,
                                    names = FALSE, type = 8)
  if (!is.finite(high_reference) || high_reference <= 0) high_reference <- scale_value
  prominence <- pmin(1, positive_change / max(high_reference, 0.25 * scale_value))
  if (is.finite(start)) prominence[x < start] <- 0

  fade_start <- if (is.finite(decline_onset)) decline_onset else NA_real_
  fade_end <- if (is.finite(decline_onset)) max(x) else NA_real_
  fade <- if (is.finite(decline_onset)) {
    ifelse(x <= fade_start, 1, pmax(0.12, 1 - 0.88 * (x - fade_start) / max(fade_end - fade_start, 1e-6)))
  } else {
    rep(1, length(x))
  }
  activity <- pmin(1, pmax(0, prominence * fade))
  alpha <- pmin(0.78, 0.04 + 0.68 * activity)
  half_height <- 0.06 + 0.20 * activity
  data.table(pseudotime = x, fitted = y, baseline = baseline,
             programme_prominence = activity, alpha = alpha,
             half_height = half_height, ymin = NA_real_, ymax = NA_real_,
             fade = fade, fade_start = fade_start, fade_end = fade_end)
}

classify_precedence <- function(probability) {
  if (!length(probability) || !is.finite(probability)) return("Not available")
  if (probability >= 0.80) return("Supported")
  if (probability >= 0.60) return("Partial")
  if (probability >= 0.40) return("Not resolved")
  "Opposite"
}

classify_precedence_evidence <- function(probability, valid_fraction, tie_fraction,
                                         delta_q025, delta_q975,
                                         onset_finite_fraction = NA_real_,
                                         tolerance = figure5_temporal_parameters$precedence_tolerance) {
  if (!is.finite(probability) || !is.finite(valid_fraction) || valid_fraction <= 0) return("Not available")
  if (valid_fraction < 0.80 || (is.finite(onset_finite_fraction) && onset_finite_fraction < 0.80)) return("Unstable")
  if (is.finite(tie_fraction) && tie_fraction >= 0.25) return("Not resolved")
  if (is.finite(delta_q025) && is.finite(delta_q975) &&
      delta_q025 <= -tolerance && delta_q975 >= tolerance) return("Not resolved")
  classify_precedence(probability)
}

patient_eligibility <- function(n_cells, n_bins, n_states, patient_meta_eligible) {
  as.logical(n_cells >= 50 & n_bins >= 3 & n_states >= 2 & patient_meta_eligible)
}

figure5_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

figure5_paths <- function(root = figure5_project_root(), create = TRUE) {
  paths <- list(
    root = root,
    scripts = file.path(root, "scripts"),
    metadata = file.path(root, "metadata", "driver", "figure5_temporal_positioning"),
    processed = file.path(root, "data", "processed", "driver", "figure5_temporal_positioning"),
    figures = file.path(root, "figures", "driver"),
    preview = file.path(root, "figures", "driver", "figure5_temporal_positioning_preview"),
    reports = file.path(root, "reports")
  )
  if (create) {
    invisible(lapply(paths[c("metadata", "processed", "preview", "reports")], dir.create,
                     recursive = TRUE, showWarnings = FALSE))
  }
  paths
}

figure5_onset_fix_paths <- function(root = figure5_project_root(), create = TRUE) {
  paths <- list(
    root = root,
    scripts = file.path(root, "scripts"),
    metadata = file.path(root, "metadata", "driver", "figure5_temporal_positioning_onset_fix"),
    processed = file.path(root, "data", "processed", "driver", "figure5_temporal_positioning_onset_fix"),
    figures = file.path(root, "figures", "driver", "figure5_temporal_positioning_onset_fix"),
    preview = file.path(root, "figures", "driver", "figure5_temporal_positioning_onset_fix_preview"),
    reports = file.path(root, "reports")
  )
  if (create) {
    invisible(lapply(paths[c("metadata", "processed", "figures", "preview", "reports")], dir.create,
                     recursive = TRUE, showWarnings = FALSE))
  }
  paths
}

# Independent output namespace for the refactored six-panel main figure.  The
# original Figure 5 correction namespace remains read-only for this workflow;
# every new panel, source-data table and preview is written below this path.
figure5_six_panel_paths <- function(root = figure5_project_root(), create = TRUE) {
  paths <- list(
    root = root,
    scripts = file.path(root, "scripts"),
    metadata = file.path(root, "metadata", "driver", "figure5_temporal_positioning_six_panel"),
    processed = file.path(root, "data", "processed", "driver", "figure5_temporal_positioning_six_panel"),
    figures = file.path(root, "figures", "driver", "figure5_temporal_positioning_six_panel"),
    preview = file.path(root, "figures", "driver", "figure5_temporal_positioning_six_panel_preview"),
    extended_metadata = file.path(root, "metadata", "driver", "figure5_temporal_positioning_six_panel", "extended_data"),
    extended_processed = file.path(root, "data", "processed", "driver", "figure5_temporal_positioning_six_panel", "extended_data"),
    extended_figures = file.path(root, "figures", "extended_data", "figure5_temporal_positioning_six_panel"),
    reports = file.path(root, "reports")
  )
  if (create) {
    invisible(lapply(paths[c("metadata", "processed", "figures", "preview",
                             "extended_metadata", "extended_processed", "extended_figures", "reports")],
                     dir.create, recursive = TRUE, showWarnings = FALSE))
  }
  paths
}

figure5_cli <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      value <- sub("^--[^=]+=", "", token)
    } else if (grepl("^--", token) && i < length(args)) {
      key <- sub("^--", "", token)
      value <- args[[i + 1L]]
      i <- i + 1L
    } else {
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", key)
    if (key %in% names(defaults)) {
      template <- defaults[[key]]
      out[[key]] <- if (is.integer(template)) as.integer(value) else if (is.numeric(template)) as.numeric(value) else if (is.logical(template)) tolower(value) %in% c("1", "true", "yes") else value
    } else {
      out[[key]] <- value
    }
    i <- i + 1L
  }
  out
}

figure5_write_json <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(object, path, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
  invisible(path)
}

figure5_write_tsv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  dt <- data.table::as.data.table(object)
  # Keep multiline figure labels readable while preserving a valid one-record-
  # per-line TSV source-data contract.
  for (column in names(dt)) {
    if (is.character(dt[[column]])) {
      dt[[column]] <- gsub("\\r?\\n", "\\\\n", dt[[column]])
    }
  }
  data.table::fwrite(dt, path, sep = "\t", quote = FALSE, na = "")
  invisible(path)
}

figure5_package_versions <- function(packages) {
  setNames(lapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else "Not available"
  }), packages)
}

figure5_scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) <= .Machine$double.eps) return(rep(NA_real_, length(x)))
  (x - r[[1]]) / diff(r)
}

figure5_axis_labels <- c(
  identity_loss = "Identity loss",
  stress_transition = "Stress transition",
  sox4_stabilization = "SOX4 stabilization"
)

figure5_axis_score_columns <- c(
  identity_loss = "identity_loss_score",
  stress_transition = "stress_transition_score",
  sox4_stabilization = "sox4_stabilization_score"
)

select_figure5_primary_pseudobulk <- function(pseudobulk) {
  dt <- as.data.table(pseudobulk)
  patient_candidate <- dt[
    aggregation_unit == "patient" & method == "main/consensus pseudotime" & eligible_patient == TRUE
  ]
  patient_coverage <- if (nrow(patient_candidate)) range(patient_candidate$pseudotime, na.rm = TRUE) else c(NA_real_, NA_real_)
  patient_coverage_ok <- nrow(patient_candidate) >= 8L && uniqueN(patient_candidate$patient_id) >= 2L &&
    all(is.finite(patient_coverage)) && patient_coverage[[1L]] <= figure5_temporal_parameters$baseline_end &&
    patient_coverage[[2L]] >= 0.90
  if (patient_coverage_ok) {
    return(list(
      data = patient_candidate,
      unit = "patient",
      analysis_unit_note = "eligible inferred-patient pseudobulk",
      fallback_reason = "",
      coverage = patient_coverage
    ))
  }

  sample_candidate <- copy(dt[aggregation_unit == "sample" & method == "main/consensus pseudotime"])
  if ("sample_id" %in% names(sample_candidate)) sample_candidate[, patient_id := sample_id]
  sample_coverage <- if (nrow(sample_candidate)) range(sample_candidate$pseudotime, na.rm = TRUE) else c(NA_real_, NA_real_)
  list(
    data = sample_candidate,
    unit = "sample",
    analysis_unit_note = "sample-token-level pseudobulk selected by the pre-specified coverage fallback",
    fallback_reason = sprintf("eligible patient-token pseudobulk coverage %.3f-%.3f failed the required <=0.10/>=0.90 rule",
                              patient_coverage[[1L]], patient_coverage[[2L]]),
    coverage = sample_coverage
  )
}

fit_landmark_gam <- function(data, score_col, k = 5L, adjusted = FALSE, use_random_effect = TRUE) {
  work <- as.data.table(data)[is.finite(pseudotime) & is.finite(get(score_col))]
  if (nrow(work) < 8L || uniqueN(work$pseudotime) < 5L) return(NULL)
  work[, dataset_id := factor(dataset_id)]
  work[, patient_id := factor(patient_id)]
  k_use <- max(3L, min(as.integer(k), uniqueN(work$pseudotime) - 1L))
  terms <- c(sprintf("s(pseudotime, k=%d, bs='cr')", k_use))
  if (uniqueN(work$dataset_id) > 1L) terms <- c(terms, "dataset_id")
  if (adjusted && "proliferation_score" %in% names(work)) terms <- c(terms, "proliferation_score")
  if (use_random_effect && uniqueN(work$patient_id) > 2L) terms <- c(terms, "s(patient_id, bs='re')")
  formula <- stats::as.formula(paste(score_col, "~", paste(terms, collapse = " + ")))
  weights <- if ("n_cells" %in% names(work)) sqrt(pmax(work$n_cells, 1)) else rep(1, nrow(work))
  tryCatch(mgcv::gam(formula, data = work, method = "REML", weights = weights), error = function(e) NULL)
}

predict_landmark_gam <- function(model, data, grid = seq(0, 1, length.out = 201L)) {
  if (is.null(model)) return(data.table(pseudotime = grid, fit = NA_real_, se = NA_real_))
  model_data <- as.data.table(data)
  newdata <- data.table(pseudotime = grid)
  if ("dataset_id" %in% names(model$model)) {
    newdata[, dataset_id := factor(levels(model$model$dataset_id)[[1]], levels = levels(model$model$dataset_id))]
  }
  if ("patient_id" %in% names(model$model)) {
    newdata[, patient_id := factor(levels(model$model$patient_id)[[1]], levels = levels(model$model$patient_id))]
  }
  if ("proliferation_score" %in% names(model$model)) {
    newdata[, proliferation_score := stats::median(model_data$proliferation_score, na.rm = TRUE)]
  }
  smooth_labels <- if (length(model$smooth)) vapply(model$smooth, function(x) x$label, character(1)) else character()
  exclude <- smooth_labels[grepl("^s\\(patient_id", smooth_labels)]
  pred <- tryCatch(stats::predict(model, newdata = newdata, se.fit = TRUE, exclude = exclude), error = function(e) NULL)
  if (is.null(pred)) return(data.table(pseudotime = grid, fit = NA_real_, se = NA_real_))
  data.table(pseudotime = grid, fit = as.numeric(pred$fit), se = as.numeric(pred$se.fit))
}

landmarks_from_table <- function(data, score_col, axis, k = 5L, adjusted = FALSE, use_random_effect = TRUE) {
  work <- as.data.table(data)[is.finite(pseudotime) & is.finite(get(score_col))]
  observed_range <- if (nrow(work)) range(work$pseudotime, na.rm = TRUE) else c(NA_real_, NA_real_)
  model <- fit_landmark_gam(work, score_col, k = k, adjusted = adjusted, use_random_effect = use_random_effect)
  pred <- predict_landmark_gam(model, work, grid = seq(0, 1, length.out = figure5_temporal_parameters$grid_n))
  raw_baseline <- work[pseudotime >= figure5_temporal_parameters$baseline_start - 1e-8 &
                         pseudotime <= figure5_temporal_parameters$baseline_end + 1e-8, get(score_col)]
  baseline <- resolve_temporal_baseline_scale(raw_baseline, pred$fit)
  landmarks <- estimate_temporal_landmarks(
    pred$pseudotime,
    pred$fit,
    axis = axis,
    baseline_scale = baseline$scale,
    baseline_scale_source = baseline$source,
    observed_range = observed_range
  )
  pred[, `:=`(baseline_scale = baseline$scale, baseline_scale_source = baseline$source,
              observed_coverage_min = observed_range[[1L]], observed_coverage_max = observed_range[[2L]])]
  list(model = model, predictions = pred, landmarks = landmarks,
       diagnostics = attr(landmarks, "diagnostics"), baseline = baseline)
}

precedence_pairs <- data.table(
  comparison = c("A before B", "B before C", "A before C"),
  upstream_axis = c("identity_loss", "stress_transition", "identity_loss"),
  downstream_axis = c("stress_transition", "sox4_stabilization", "sox4_stabilization")
)

figure5_result_text <- function(precedence) {
  p <- setNames(precedence$probability, precedence$comparison)
  ab <- unname(p[["A before B"]]); bc <- unname(p[["B before C"]]); ac <- unname(p[["A before C"]])
  if (all(is.finite(c(ab, bc, ac))) && all(c(ab, bc, ac) >= 0.80)) {
    return("Pseudotemporal analyses support an ordered progression from hepatocyte identity loss through stress-transition activation to SOX4-associated malignant-state stabilization.")
  }
  if (is.finite(ac) && is.finite(bc) && ac >= 0.80 && bc >= 0.80 && (!is.finite(ab) || ab < 0.80)) {
    return("Pseudotemporal analyses support earlier hepatocyte identity loss and stress-transition activity relative to SOX4-associated malignant-state stabilization, whereas the ordering between identity loss and stress activation remains unresolved.")
  }
  if (is.finite(bc) && bc >= 0.80) {
    return("AP-1/CEBPB/EGR1 stress-transition activity tended to precede SOX4-associated malignant-state stabilization, while the relative timing of HNF4A/PPARA identity loss was not consistently resolved.")
  }
  "The three regulatory programmes occupied overlapping regions of the hepatocyte-state continuum, but their strict temporal ordering could not be resolved."
}
