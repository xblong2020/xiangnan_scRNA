#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_onset_fix_paths(root, create = FALSE)
boot <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
keep <- c("onset_time", "t10", "t50", "maximum_slope_time", "peak_time")
plot_data <- boot[landmark %chin% keep, {
  finite <- time[is.finite(time)]
  .(median = if (length(finite)) median(finite) else NA_real_,
    lower = if (length(finite)) quantile(finite, 0.025, type = 8) else NA_real_,
    upper = if (length(finite)) quantile(finite, 0.975, type = 8) else NA_real_,
    n_bootstrap = .N, n_finite = length(finite), finite_fraction = length(finite) / .N)
}, by = .(axis, landmark)]
plot_data[, stability := ifelse(landmark %chin% c("onset_time", "t10") & finite_fraction < 0.80, "Unstable", "Stable")]
plot_data[, landmark_label := factor(landmark, levels = keep,
                                     labels = c("Corrected onset", "t10", "t50 (core)", "Maximum slope (core)", "Peak"))]
axis_levels <- rev(names(figure5_axis_score_columns))
plot_data[, y_base := match(axis, axis_levels)]
offsets <- setNames(seq(-0.24, 0.24, length.out = length(levels(plot_data$landmark_label))), levels(plot_data$landmark_label))
plot_data[, y := y_base + offsets[as.character(landmark_label)]]

p <- ggplot(plot_data, aes(colour = axis, shape = landmark_label)) +
  geom_segment(data = plot_data[is.finite(median)],
               aes(x = lower, xend = upper, y = y, yend = y, linetype = stability), linewidth = 0.55) +
  geom_point(data = plot_data[is.finite(median)], aes(x = median, y = y), size = 2.2, fill = "white") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  scale_shape_manual(values = c("Corrected onset" = 16, "t10" = 1, "t50 (core)" = 18,
                                "Maximum slope (core)" = 17, "Peak" = 15)) +
  scale_linetype_manual(values = c(Stable = "solid", Unstable = "dashed"), guide = "none") +
  scale_y_continuous(breaks = seq_along(axis_levels), labels = figure5_axis_labels[axis_levels],
                     expand = expansion(add = 0.45)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "D  Corrected bootstrap temporal landmarks",
       subtitle = "Corrected onset, t10, t50, maximum slope and peak; peak is retained in D only",
       x = "Oriented consensus pseudotime", y = NULL, shape = "Landmark") +
  theme_figure5() + theme(legend.position = "bottom")

out_dir <- file.path(paths$figures, "panel_D_landmarks")
outputs <- export_figure5_plot(p, file.path(out_dir, "figure5_six_panel_D_landmarks"), 5.9, 3.5)
figure5_write_tsv(plot_data, file.path(paths$metadata, "figure5_six_panel_D_temporal_landmarks.tsv"))
figure5_write_json(list(
  panel = "5D",
  title = "Corrected bootstrap temporal landmarks",
  interval = "2.5th-97.5th bootstrap percentiles",
  retained_landmarks = c("corrected onset", "t10", "t50", "maximum slope", "peak"),
  core_markers = c("t50", "maximum_slope_time"),
  bootstrap_unit = unique(boot$bootstrap_unit),
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_D_report.json"))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5_panel_D.rds"))

