#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))

expect_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

expect_na <- function(value, message) {
  if (!(length(value) == 1L && is.na(value))) stop(message, call. = FALSE)
}

fixed <- list(
  baseline_end = 0.10,
  search_start = 0.10,
  slope_search_end = 0.95,
  grid_n = 201L,
  onset_min_run = 5L,
  crossing_min_run = 3L,
  min_effect_fraction = 0.05,
  min_effect_baseline_sd = 0.25,
  derivative_fraction = 0.10,
  precedence_tolerance = 0.005
)

x <- seq(0, 1, length.out = fixed$grid_n)

# 1. Flat curve.
flat <- estimate_temporal_landmarks(x, rep(2, length(x)), baseline_scale = 0.1)
expect_true(identical(names(flat), c("onset_time", "t10", "t50", "maximum_slope_time",
                                     "extremum_time", "peak_time", "plateau_time", "decline_onset")),
            "Landmark API must retain exactly the eight contract fields")
for (field in c("onset_time", "t10", "t50", "maximum_slope_time")) {
  expect_na(flat[[field]], paste("Flat curve must not estimate", field))
}

# 2. Delayed sigmoid centred at 0.35.
sigmoid <- 1 / (1 + exp(-(x - 0.35) / 0.035))
sig <- estimate_temporal_landmarks(x, sigmoid, baseline_scale = 0.02)
expect_true(is.finite(sig$onset_time) && sig$onset_time > 0.10 && sig$onset_time < 0.40,
            "Delayed sigmoid onset must lie in the post-baseline rise")
expect_true(is.finite(sig$maximum_slope_time) && abs(sig$maximum_slope_time - 0.35) <= 0.03,
            "Delayed sigmoid maximum slope must recover the true rise centre")

# 3. Early negative dip followed by a late rise.
early_dip <- -0.35 * exp(-((x - 0.12) / 0.035)^2) + 1 / (1 + exp(-(x - 0.48) / 0.04))
dip <- estimate_temporal_landmarks(x, early_dip, baseline_scale = 0.03)
expect_true(is.finite(dip$onset_time) && dip$onset_time > 0.25,
            "Early reverse-direction dip must not be called activation onset")

# 4. Positive slope confined to the baseline boundary window.
boundary_only <- pmin(x, 0.10)
boundary <- estimate_temporal_landmarks(x, boundary_only, baseline_scale = 0.01)
expect_na(boundary$onset_time, "Baseline-window slope must not create onset")
expect_na(boundary$maximum_slope_time, "Baseline-window slope must not create maximum slope")

# 5. Transient rise shorter than five grid points.
short_spike <- rep(0, length(x))
short_spike[45:48] <- c(0.2, 0.6, 0.6, 0.2)
spike <- estimate_temporal_landmarks(x, short_spike, baseline_scale = 0.02)
expect_na(spike$onset_time, "A transient shorter than five grid points must not create onset")

# 6. Monotonic rise starting after 0.15.
monotonic <- pmax(0, x - 0.15)
mono <- estimate_temporal_landmarks(x, monotonic, baseline_scale = 0.02)
expect_true(is.finite(mono$onset_time) && mono$onset_time > 0.10,
            "Post-0.15 monotonic rise must have a post-baseline onset")
expect_true(is.finite(mono$t10) && is.finite(mono$t50) && mono$onset_time <= mono$t10 && mono$t10 <= mono$t50,
            "Monotonic landmarks must follow onset <= t10 <= t50")

accelerating_to_boundary <- x^3
boundary_limited <- estimate_temporal_landmarks(x, accelerating_to_boundary, baseline_scale = 0.02)
expect_na(boundary_limited$maximum_slope_time,
          "A maximum derivative attained only at the right search boundary must be unresolved")
expect_true(isTRUE(attr(boundary_limited, "diagnostics")$maximum_slope_boundary_hit),
            "Right-boundary maximum slope must be diagnosed")

# 7. Bell-shaped stress curve.
bell <- exp(-((x - 0.60) / 0.16)^2)
stress <- estimate_temporal_landmarks(x, bell, baseline_scale = 0.02)
expect_true(is.finite(stress$onset_time) && stress$onset_time < stress$peak_time,
            "Bell-shaped curve must have onset before peak")
expect_true(is.finite(stress$peak_time) && abs(stress$peak_time - 0.60) <= 0.02,
            "Bell-shaped curve peak must be recovered")
expect_true(is.finite(stress$decline_onset) && stress$decline_onset > stress$peak_time,
            "Bell-shaped curve must have a sustained decline after peak")

# 8. Tied temporal landmarks contribute 0.5.
tied <- tie_aware_precedence_probability(c(0.30, 0.40), c(0.30, 0.40), tolerance = fixed$precedence_tolerance)
expect_true(identical(tied$n_tied, 2L) && tied$probability == 0.5,
            "Tied landmarks must contribute 0.5 rather than zero")

# 9. Missing early coverage.
x_missing_early <- seq(0.20, 1, length.out = 161)
missing_early <- estimate_temporal_landmarks(x_missing_early, plogis((x_missing_early - 0.5) / 0.05), baseline_scale = 0.02)
expect_true(all(vapply(missing_early, is.na, logical(1))), "Missing early coverage must return all NA landmarks")
expect_true(identical(attr(missing_early, "diagnostics")$coverage_ok, FALSE),
            "Missing early coverage must be recorded as a coverage failure")

# 10. Missing late coverage.
x_missing_late <- seq(0, 0.85, length.out = 171)
missing_late <- estimate_temporal_landmarks(x_missing_late, plogis((x_missing_late - 0.5) / 0.05), baseline_scale = 0.02)
expect_true(all(vapply(missing_late, is.na, logical(1))), "Missing late coverage must return all NA landmarks")
expect_true(grepl("late", attr(missing_late, "diagnostics")$failure_reason),
            "Missing late coverage must not extrapolate a stable peak")

message("PASS: 10 Figure 5 temporal-landmark synthetic tests")
