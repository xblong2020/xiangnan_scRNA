#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)

boot <- as.data.table(readRDS(file.path(paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
gam_predictions <- as.data.table(readRDS(file.path(paths$processed, "figure5b_gam_predictions.rds")))[scenario == "primary"]
precedence <- fread(file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))[is_primary == TRUE]
cells <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_oriented_cell_scores.rds")))[method == "main/consensus pseudotime"]

summary_landmark <- function(axis_name, landmark_name) {
  values <- boot[axis == axis_name & landmark == landmark_name, time]
  finite <- values[is.finite(values)]
  if (!length(finite)) return(c(median = NA_real_, lower = NA_real_, upper = NA_real_, n = length(values), finite_fraction = 0))
  c(
    median = stats::median(finite),
    lower = unname(stats::quantile(finite, 0.025, type = 8)),
    upper = unname(stats::quantile(finite, 0.975, type = 8)),
    n = length(values),
    finite_fraction = length(finite) / length(values)
  )
}

is_stable_landmark <- function(summary, minimum_fraction = 0.80, maximum_width = 0.25) {
  all(is.finite(summary[c("median", "lower", "upper", "finite_fraction")])) &&
    summary[["finite_fraction"]] >= minimum_fraction &&
    summary[["upper"]] - summary[["lower"]] <= maximum_width
}

bands <- rbindlist(lapply(names(figure5_axis_score_columns), function(axis_name) {
  onset <- summary_landmark(axis_name, "onset_time")
  t10 <- summary_landmark(axis_name, "t10")
  t50 <- summary_landmark(axis_name, "t50")
  maximum_slope <- summary_landmark(axis_name, "maximum_slope_time")
  decline <- summary_landmark(axis_name, "decline_onset")
  choice <- resolve_figure5h_start(
    onset = onset[["median"]], onset_lower = onset[["lower"]], onset_upper = onset[["upper"]],
    onset_finite_fraction = onset[["finite_fraction"]],
    t10 = t10[["median"]], t10_lower = t10[["lower"]], t10_upper = t10[["upper"]],
    t10_finite_fraction = t10[["finite_fraction"]]
  )
  stable_decline <- is_stable_landmark(decline)
  prediction <- gam_predictions[axis == axis_name]
  profile <- build_figure5h_activity_profile(
    pseudotime = prediction$pseudotime,
    fitted = prediction$fit,
    start = choice$start,
    t50 = t50[["median"]],
    maximum_slope = maximum_slope[["median"]],
    decline_onset = if (stable_decline) decline[["median"]] else NA_real_,
    direction = 1
  )
  profile[, `:=`(
    axis = axis_name,
    y = match(axis_name, rev(names(figure5_axis_score_columns))),
    axis_label = figure5_axis_labels[[axis_name]],
    onset = onset[["median"]], onset_lower = onset[["lower"]], onset_upper = onset[["upper"]],
    onset_finite_fraction = onset[["finite_fraction"]], onset_stable = is_stable_landmark(onset),
    t10 = t10[["median"]], t10_lower = t10[["lower"]], t10_upper = t10[["upper"]],
    t10_finite_fraction = t10[["finite_fraction"]], t10_stable = is_stable_landmark(t10),
    t50 = t50[["median"]], t50_lower = t50[["lower"]], t50_upper = t50[["upper"]],
    maximum_slope = maximum_slope[["median"]],
    maximum_slope_lower = maximum_slope[["lower"]], maximum_slope_upper = maximum_slope[["upper"]],
    decline_onset = decline[["median"]], decline_lower = decline[["lower"]], decline_upper = decline[["upper"]],
    decline_finite_fraction = decline[["finite_fraction"]], decline_stable = stable_decline,
    boundary_start = choice$start, boundary_start_lower = choice$start_lower,
    boundary_start_upper = choice$start_upper, boundary_method = choice$method,
    boundary_status = choice$boundary_status,
    right_edge_rule = if (stable_decline) "fade after stable decline_onset" else "continued prominence to observed pseudotime end"
  )]
  profile[, `:=`(ymin = y - half_height, ymax = y + half_height,
                 alpha = pmin(0.78, pmax(0.04, alpha)))]
  profile
}), fill = TRUE)

band_summary <- unique(bands[, .(
  axis, axis_label, y,
  onset, onset_lower, onset_upper, onset_finite_fraction, onset_stable,
  t10, t10_lower, t10_upper, t10_finite_fraction, t10_stable,
  t50, t50_lower, t50_upper,
  maximum_slope, maximum_slope_lower, maximum_slope_upper,
  decline_onset, decline_lower, decline_upper, decline_finite_fraction, decline_stable,
  boundary_start, boundary_start_lower, boundary_start_upper, boundary_method, boundary_status,
  fade_start, fade_end, right_edge_rule
)])
band_summary[, unresolved_x := 0.13]

band_polygons <- rbindlist(lapply(split(bands, bands$axis), function(values) {
  if (nrow(values) < 2L) return(NULL)
  rbindlist(lapply(seq_len(nrow(values) - 1L), function(i) {
    data.table(
      axis = values$axis[[i]],
      segment_group = paste(values$axis[[i]], i, sep = "::"),
      x = c(values$pseudotime[[i]], values$pseudotime[[i + 1L]], values$pseudotime[[i + 1L]], values$pseudotime[[i]]),
      y = c(values$ymin[[i]], values$ymin[[i + 1L]], values$ymax[[i + 1L]], values$ymax[[i]]),
      alpha = mean(c(values$alpha[[i]], values$alpha[[i + 1L]]), na.rm = TRUE)
    )
  }))
}), fill = TRUE)

outline <- rbindlist(lapply(split(bands, bands$axis), function(values) rbind(
  values[, .(axis, pseudotime, outline_y = ymax + 0.018, group = paste(axis, "upper", sep = "::"))],
  values[, .(axis, pseudotime, outline_y = ymin - 0.018, group = paste(axis, "lower", sep = "::"))]
)), fill = TRUE)

cells[, state_group := fcase(
  trajectory_role == "normal_reference", "Normal/reference",
  grepl("stress", trajectory_role, ignore.case = TRUE), "Stressed/injured",
  grepl("regener", trajectory_role, ignore.case = TRUE), "Regenerative/progenitor",
  grepl("prolifer", trajectory_role, ignore.case = TRUE), "Proliferating candidate",
  grepl("malignant", trajectory_role, ignore.case = TRUE), "Malignant/malignant-like",
  default = "Stressed/injured"
)]
state_levels <- c("Normal/reference", "Stressed/injured", "Regenerative/progenitor", "Proliferating candidate", "Malignant/malignant-like")
state_med <- cells[, .(state_marker = median(pseudotime, na.rm = TRUE), n_cells = .N), by = state_group]
state_med <- merge(data.table(state_group = state_levels, fallback = seq(0.08, 0.92, length.out = 5)), state_med, by = "state_group", all.x = TRUE)
state_med[, marker_x := fifelse(is.finite(state_marker), state_marker, fallback)]
state_med[, state_order := match(state_group, state_levels)]
setorder(state_med, state_order)
state_med[, label_y := ifelse(state_order %% 2L == 0L, 3.82, 3.55)]
state_med[, state_label := fcase(
  state_group == "Normal/reference", "Normal/\nreference",
  state_group == "Stressed/injured", "Stressed/\ninjured",
  state_group == "Regenerative/progenitor", "Regenerative/\nprogenitor",
  state_group == "Proliferating candidate", "Proliferating\ncandidate",
  default = "Malignant/\nmalignant-like"
)]

overall_support <- if (all(precedence$evidence_grade == "Supported")) "Supported" else if (any(precedence$evidence_grade %chin% c("Partial", "Unstable", "Not resolved", "Opposite"))) "partial or unresolved" else "not resolved"
caption_text <- "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events."

p <- ggplot() +
  geom_segment(data = state_med, aes(x = marker_x, xend = marker_x, y = 3.25, yend = 3.42),
               colour = lancet_palette[8], linewidth = 0.35, linetype = "dotted") +
  geom_text(data = state_med, aes(x = marker_x, y = label_y, label = state_label),
            size = 1.9, lineheight = 0.88, colour = lancet_palette[9]) +
  geom_polygon(data = band_polygons, aes(x = x, y = y, fill = axis, alpha = alpha, group = segment_group), colour = NA) +
  geom_line(data = outline, aes(x = pseudotime, y = outline_y, group = group),
            colour = lancet_palette[8], linewidth = 0.30, linetype = "dashed", alpha = 0.55) +
  geom_segment(data = band_summary[is.finite(boundary_start)],
               aes(x = boundary_start_lower, xend = boundary_start_upper, y = y, yend = y),
               colour = lancet_palette[8], linewidth = 0.55, linetype = "dashed") +
  geom_point(data = band_summary[is.finite(boundary_start)], aes(x = boundary_start, y = y, colour = axis),
             size = 2.0, shape = 21, fill = "white", stroke = 0.55) +
  geom_segment(data = band_summary[is.finite(t50)], aes(x = t50_lower, xend = t50_upper, y = y - 0.05, yend = y - 0.05),
               colour = lancet_palette[8], linewidth = 0.35, linetype = "dashed") +
  geom_point(data = band_summary[is.finite(t50)], aes(x = t50, y = y - 0.05, colour = axis), size = 2.0, shape = 16) +
  geom_segment(data = band_summary[is.finite(maximum_slope)],
               aes(x = maximum_slope_lower, xend = maximum_slope_upper, y = y + 0.05, yend = y + 0.05),
               colour = lancet_palette[8], linewidth = 0.35, linetype = "dashed") +
  geom_point(data = band_summary[is.finite(maximum_slope)], aes(x = maximum_slope, y = y + 0.05, colour = axis), size = 2.2, shape = 17) +
  geom_text(data = band_summary, aes(x = 0.01, y = y + 0.30, label = axis_label, colour = axis),
            hjust = 0, fontface = "bold", size = 2.45) +
  geom_text(data = band_summary[boundary_status == "boundary unresolved"],
            aes(x = unresolved_x, y = y, label = "boundary unresolved"),
            hjust = 0, size = 2.0, colour = lancet_palette[8]) +
  annotate("text", x = 0.5, y = 0.42, label = caption_text,
           fontface = "italic", size = 2.10, colour = lancet_palette[9]) +
  annotate("text", x = 0.5, y = 0.16, label = "Pseudotemporal phases overlap and do not establish a strict causal cascade.",
           fontface = "italic", size = 2.30, colour = lancet_palette[9]) +
  annotate("text", x = 0.99, y = 0.62, label = overall_support, hjust = 1, size = 2.20, colour = lancet_palette[8]) +
  scale_fill_manual(values = axis_palette, guide = "none") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 4.05), breaks = NULL) +
  labs(
    title = "H  Corrected overlapping regulatory-activity model",
    subtitle = "Open circle: resolved start; circle: t50; triangle: maximum slope; bands fade only after stable decline",
    x = "Relative hepatocyte-state progression", y = NULL
  ) +
  theme_figure5() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 7.2, colour = lancet_palette[8]))

boundaries <- rbindlist(list(
  band_summary[, .(
    record_type = "axis_activity_band", axis, label = axis_label,
    onset, onset_lower, onset_upper, onset_finite_fraction, onset_stable,
    t10, t10_lower, t10_upper, t10_finite_fraction, t10_stable,
    t50, t50_lower, t50_upper,
    maximum_slope, maximum_slope_lower, maximum_slope_upper,
    decline_onset, decline_lower, decline_upper, decline_finite_fraction, decline_stable,
    boundary_start, boundary_start_lower, boundary_start_upper, boundary_method, boundary_status,
    fade_start, fade_end, right_edge_rule, programme_end_estimated = FALSE,
    source = "1,000 patient/sample-stratified corrected REML bootstrap replicates"
  )],
  state_med[, .(
    record_type = "background_state_marker", axis = "", label = state_group,
    onset = NA_real_, onset_lower = NA_real_, onset_upper = NA_real_, onset_finite_fraction = NA_real_, onset_stable = FALSE,
    t10 = NA_real_, t10_lower = NA_real_, t10_upper = NA_real_, t10_finite_fraction = NA_real_, t10_stable = FALSE,
    t50 = NA_real_, t50_lower = NA_real_, t50_upper = NA_real_,
    maximum_slope = NA_real_, maximum_slope_lower = NA_real_, maximum_slope_upper = NA_real_,
    decline_onset = NA_real_, decline_lower = NA_real_, decline_upper = NA_real_, decline_finite_fraction = NA_real_, decline_stable = FALSE,
    boundary_start = marker_x, boundary_start_lower = marker_x, boundary_start_upper = marker_x,
    boundary_method = "observed state median marker", boundary_status = "marker only",
    fade_start = NA_real_, fade_end = NA_real_, right_edge_rule = "not applicable", programme_end_estimated = FALSE,
    source = ifelse(is.finite(state_marker), "observed median", "fallback marker")
  )]
), fill = TRUE)
figure5_write_tsv(boundaries, file.path(paths$metadata, "figure5h_phase_boundaries.tsv"))
figure5_write_tsv(band_summary, file.path(paths$metadata, "figure5h_activity_band_summary.tsv"))

out_dir <- file.path(paths$figures, "figure5h_overlapping_phase_model")
export_figure5_plot(p, file.path(out_dir, "figure5h_overlapping_phase_model"), 7.2, 4.05)
figure5_write_json(list(
  formal_band_left_boundary = "corrected onset when finite fraction >=0.80 and CI width <=0.25; otherwise stable t10; otherwise boundary unresolved",
  uses_old_onset_for_formal_band = FALSE,
  corrected_onset_retained_in_source_data = TRUE,
  right_edge = "fade only after stable decline_onset; otherwise activity prominence continues to pseudotime end",
  no_discrete_stage_rectangles = TRUE,
  unresolved_axes = band_summary[boundary_status == "boundary unresolved", axis],
  caption = caption_text,
  support = overall_support,
  causal_claim = FALSE,
  precedence = split(precedence, seq_len(nrow(precedence)))
), file.path(paths$metadata, "figure5h_overlapping_phase_model_report.json"))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5h.rds"))
