#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(patchwork))
paths <- figure5_onset_fix_paths(root)
pred <- as.data.table(readRDS(file.path(paths$processed, "figure5b_gam_predictions.rds")))
pred[, axis_label := factor(figure5_axis_labels[axis], levels = figure5_axis_labels)]

main <- pred[scenario == "primary"]
p_main <- ggplot(main, aes(pseudotime, fit, colour = axis, fill = axis)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~axis_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = axis_palette, guide = "none") + scale_fill_manual(values = axis_palette, guide = "none") +
  labs(title = "B  Coverage-corrected regulatory axes along unified pseudotime", x = "Normalized pseudotime", y = "Scaled programme score") +
  theme_figure5()

selected <- c("primary", "adjusted_proliferation", "no_high_proliferation", "no_generic_stress", "regulon_auc", "tf_expression")
overlay <- pred[scenario %chin% selected]
p_overlay <- ggplot(overlay, aes(pseudotime, fit, colour = axis, linetype = scenario, group = interaction(axis, scenario))) +
  geom_line(linewidth = 0.55, alpha = 0.9) + facet_wrap(~axis_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = axis_palette) +
  labs(title = "Figure 5B sensitivity overlay", x = "Normalized pseudotime", y = "Scaled programme score", colour = "Axis", linetype = "Analysis") +
  theme_figure5() + theme(legend.position = "bottom")

make_single <- function(scenario_name, title) {
  ggplot(pred[scenario == scenario_name], aes(pseudotime, fit, colour = axis, fill = axis)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) + geom_line(linewidth = 0.75) +
    facet_wrap(~axis_label, ncol = 1, scales = "free_y") + scale_colour_manual(values = axis_palette, guide = "none") +
    scale_fill_manual(values = axis_palette, guide = "none") + labs(title = title, x = "Normalized pseudotime", y = "Scaled programme score") + theme_figure5()
}
p_adjusted <- make_single("adjusted_proliferation", "Figure 5B adjusted for proliferation")
p_cnv <- make_single("cnv_strict", "Figure 5B CNV-strict sensitivity")

out_dir <- file.path(paths$figures, "figure5b_three_axis_gam")
export_figure5_plot(p_main, file.path(out_dir, "figure5b_three_axis_gam_main"), 3.7, 6.4)
export_figure5_plot(p_overlay, file.path(out_dir, "figure5b_three_axis_gam_overlay_sensitivity"), 5.4, 6.4)
export_figure5_plot(p_adjusted, file.path(out_dir, "figure5b_three_axis_gam_adjusted"), 3.7, 6.4)
export_figure5_plot(p_cnv, file.path(out_dir, "figure5b_three_axis_gam_cnv_strict"), 3.7, 6.4)
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p_main, file.path(paths$processed, "plot_objects", "figure5b.rds"))
