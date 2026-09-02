#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_figure5_temporal_positioning.R"
root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"), local = FALSE)

expect_true <- function(value, label) {
  if (!isTRUE(value)) stop(label, call. = FALSE)
}

expect_equal <- function(actual, expected, tolerance = 1e-8, label = "values differ") {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tolerance, na.rm = TRUE)) {
    stop(sprintf("%s: actual=%s expected=%s", label, paste(actual, collapse = ","), paste(expected, collapse = ",")), call. = FALSE)
  }
}

# Explicit sample-name rules must merge paired material from the same patient.
mapped <- derive_figure5_patient_id(
  dataset_id = c("GSE149614", "GSE149614", "GSE185477", "GSE185477", "GSE202379"),
  cnv_sample = c("HCC03T", "HCC03N", "C41_TST", "C41_NST", "GSE202379_SeuratObject_AllCells_counts_fixed")
)
expect_true(mapped$patient_id[[1]] == mapped$patient_id[[2]], "HCC03 paired samples were not merged")
expect_true(mapped$patient_id[[3]] == mapped$patient_id[[4]], "C41 paired samples were not merged")
expect_true(mapped$patient_id_source[[5]] == "aggregate_sample_proxy", "aggregate object must remain explicitly flagged")
expect_true(!mapped$patient_meta_eligible[[5]], "aggregate object must be excluded from patient meta-analysis")

# Dataset-wise robust scaling must center non-constant groups and keep constants finite.
scaled <- robust_z_by_group(c(1, 2, 3, 8, 8, 8), c("a", "a", "a", "b", "b", "b"))
expect_equal(stats::median(scaled[1:3]), 0, label = "robust z-score was not centered")
expect_equal(scaled[4:6], c(0, 0, 0), label = "constant group was not mapped to zero")

# A reversed trajectory must be detected from endpoint evidence.
pt <- seq(0, 1, length.out = 101)
orientation <- decide_pseudotime_orientation(
  pseudotime = rev(pt),
  mature_identity_score = 1 - pt,
  cnv_score = pt,
  malignant_fate = pt,
  malignant_indicator = as.numeric(pt > 0.7)
)
expect_true(isTRUE(orientation$flipped), "reversed pseudotime was not flipped")
expect_true(orientation$oriented_score > 0, "oriented endpoint score must be positive")

# Landmark estimates must follow the known synthetic ordering.
grid <- seq(0, 1, length.out = 401)
axis_a <- 1 / (1 + exp(-24 * (grid - 0.28)))
axis_b <- exp(-((grid - 0.52) / 0.16)^2)
axis_c <- 1 / (1 + exp(-22 * (grid - 0.72)))
la <- estimate_temporal_landmarks(grid, axis_a, axis = "identity_loss")
lb <- estimate_temporal_landmarks(grid, axis_b, axis = "stress_transition")
lc <- estimate_temporal_landmarks(grid, axis_c, axis = "sox4_stabilization")
expect_true(abs(la$maximum_slope_time - 0.28) < 0.03, "Axis A maximum slope is inaccurate")
expect_true(abs(lb$peak_time - 0.52) < 0.03, "Axis B peak is inaccurate")
expect_true(abs(lc$t50 - 0.72) < 0.04, "Axis C t50 is inaccurate")
expect_true(la$maximum_slope_time < lb$maximum_slope_time && lb$maximum_slope_time < lc$maximum_slope_time,
            "synthetic temporal ordering was not recovered")

# Figure 5H activation boundaries must reject early reverse-direction movement.
activation_grid <- seq(0, 1, length.out = 101)
activation_curve <- ifelse(activation_grid < 0.25,
                           0.35 - 0.35 * activation_grid / 0.25,
                           0.35 + 1.25 * (activation_grid - 0.25))
activation_derivative <- c(diff(activation_curve) / diff(activation_grid), NA_real_)
activation_lower <- activation_derivative - 0.15
activation_upper <- activation_derivative + 0.15
directional <- find_directional_activation(
  pseudotime = activation_grid,
  fitted = activation_curve,
  derivative = activation_derivative,
  derivative_lower = activation_lower,
  derivative_upper = activation_upper,
  direction = 1,
  min_change_sd = 0.25,
  min_run = 5L
)
expect_true(isTRUE(directional$found), "directional activation was not detected")
expect_true(directional$time >= 0.25, "early reverse decline was incorrectly treated as activation")

t10_choice <- resolve_figure5h_start(
  t10 = 0.10,
  t10_lower = 0.08,
  t10_upper = 0.13,
  t10_finite_fraction = 1,
  directional = directional,
  max_t10_ci_width = 0.25
)
expect_true(t10_choice$method == "directional_activation_fallback", "invalid directional t10 was not replaced")
expect_true(t10_choice$start >= directional$time, "formal H boundary precedes directional activation")

stable_t10 <- resolve_figure5h_start(
  t10 = 0.12,
  t10_lower = 0.10,
  t10_upper = 0.15,
  t10_finite_fraction = 1,
  directional = list(found = TRUE, time = 0.08),
  max_t10_ci_width = 0.25
)
expect_true(stable_t10$method == "bootstrap_t10", "stable t10 was not preferred")

unresolved <- resolve_figure5h_start(
  t10 = NA_real_,
  t10_lower = NA_real_,
  t10_upper = NA_real_,
  t10_finite_fraction = 0,
  directional = list(found = FALSE, time = NA_real_),
  max_t10_ci_width = 0.25
)
expect_true(unresolved$boundary_status == "boundary unresolved", "unstable boundary was not labelled")

profile <- build_figure5h_activity_profile(
  pseudotime = activation_grid,
  fitted = activation_curve,
  start = t10_choice$start,
  t50 = 0.55,
  maximum_slope = 0.55,
  direction = 1
)
expect_true(all(profile$alpha >= 0 & profile$alpha <= 1), "activity alpha is outside [0,1]")
expect_true(profile$alpha[which.min(abs(profile$pseudotime - 0.05))] < profile$alpha[which.min(abs(profile$pseudotime - 0.55))],
            "activity band is not more opaque at high activity")
expect_true(all(profile$fade == 1), "activity band was forced to fade without a stable decline")
declining_profile <- build_figure5h_activity_profile(
  pseudotime = activation_grid,
  fitted = activation_curve,
  start = t10_choice$start,
  t50 = 0.55,
  maximum_slope = 0.55,
  decline_onset = 0.70,
  direction = 1
)
expect_true(declining_profile$fade[[nrow(declining_profile)]] < 1,
            "stable decline did not produce a gradual right-edge fade")

# Evidence categories are fixed by the Figure 5 contract.
expect_true(classify_precedence(0.85) == "Supported", "Supported threshold is wrong")
expect_true(classify_precedence(0.70) == "Partial", "Partial threshold is wrong")
expect_true(classify_precedence(0.50) == "Not resolved", "Not resolved threshold is wrong")
expect_true(classify_precedence(0.30) == "Opposite", "Opposite threshold is wrong")
expect_true(classify_precedence(NA_real_) == "Not available", "NA threshold is wrong")

# Eligibility requires cells, bins, states, and a real patient identifier.
elig <- patient_eligibility(
  n_cells = c(50, 49, 100, 100),
  n_bins = c(3, 5, 2, 5),
  n_states = c(2, 3, 3, 1),
  patient_meta_eligible = c(TRUE, TRUE, TRUE, TRUE)
)
expect_true(identical(elig, c(TRUE, FALSE, FALSE, FALSE)), "patient eligibility rule is wrong")

cat("Figure 5 temporal core tests passed\n")
