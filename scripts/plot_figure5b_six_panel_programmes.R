#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_onset_fix_paths(root, create = FALSE)
pred <- as.data.table(readRDS(file.path(source_paths$processed, "figure5b_gam_predictions.rds")))
pred[, axis_label := factor(figure5_axis_labels[axis], levels = figure5_axis_labels)]

main <- pred[scenario == "primary"]
p_main <- ggplot(main, aes(pseudotime, fit, colour = axis, fill = axis)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~axis_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  scale_fill_manual(values = axis_palette, guide = "none") +
  labs(title = "B  Coverage-corrected regulatory programmes along unified pseudotime",
       subtitle = "Primary coverage-qualified sample-token pseudobulk GAM; ribbons show fitted uncertainty",
       x = "Oriented consensus pseudotime", y = "Scaled programme score") +
  theme_figure5()

selected <- c("primary", "adjusted_proliferation", "no_high_proliferation", "no_generic_stress", "regulon_auc", "tf_expression")
overlay <- pred[scenario %chin% selected]
p_overlay <- ggplot(overlay, aes(pseudotime, fit, colour = axis, linetype = scenario, group = interaction(axis, scenario))) +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  facet_wrap(~axis_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = axis_palette) +
  labs(title = "B  Coverage-corrected programme sensitivity overlay", x = "Oriented consensus pseudotime",
       y = "Scaled programme score", colour = "Programme", linetype = "Analysis") +
  theme_figure5() + theme(legend.position = "bottom")

out_dir <- file.path(paths$figures, "panel_B_programmes")
outputs <- export_figure5_plot(p_main, file.path(out_dir, "figure5_six_panel_B_programmes"), 4.2, 6.5)
export_figure5_plot(p_overlay, file.path(out_dir, "figure5_six_panel_B_programmes_sensitivity"), 5.8, 6.5)
figure5_write_tsv(pred, file.path(paths$metadata, "figure5_six_panel_B_gam_predictions.tsv.gz"))
figure5_write_json(list(
  panel = "5B",
  title = "Coverage-corrected regulatory programmes along unified pseudotime",
  analysis_unit = "Coverage-qualified patient/sample-token pseudobulk",
  primary_scenario = "primary",
  fallback_note = "Because eligible patient-token data did not span the complete pseudotemporal range, the primary bootstrap used the prespecified sample-token coverage fallback.",
  source_namespace = "figure5_temporal_positioning_onset_fix",
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_B_report.json"))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p_main, file.path(paths$processed, "plot_objects", "figure5_panel_B.rds"))

