#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_onset_fix_paths(root, create = FALSE)
dt <- as.data.table(readRDS(file.path(source_paths$processed, "figure5c_heatmap_predictions.rds")))
order_dt <- unique(dt[, .(entity, display_order, axis, entity_type)])[order(display_order)]
dt[, entity_factor := factor(entity, levels = rev(order_dt$entity))]
order_dt[, entity_factor := factor(entity, levels = rev(order_dt$entity))]

p <- ggplot(dt, aes(pseudotime, entity_factor, fill = row_z)) +
  geom_tile() +
  geom_segment(data = order_dt, aes(x = -0.035, xend = -0.012, y = entity_factor,
                                    yend = entity_factor, colour = axis),
               inherit.aes = FALSE, linewidth = 2.1) +
  scale_fill_gradientn(colours = figure5_diverging_palette, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "Row z-score") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  coord_cartesian(xlim = c(-0.04, 1), clip = "off") +
  labs(title = "C  Temporal organization of frozen TF, regulon and target-gene programmes",
       subtitle = "Frozen programme members are ordered by corrected temporal landmarks; heatmap values remain frozen",
       x = "Oriented consensus pseudotime", y = NULL) +
  theme_figure5(7) +
  theme(axis.text.y = element_text(size = 5.8), plot.margin = margin(5, 6, 5, 18))

out_dir <- file.path(paths$figures, "panel_C_programme_heatmap")
outputs <- export_figure5_plot(p, file.path(out_dir, "figure5_six_panel_C_programme_heatmap"), 5.8, max(6.0, nrow(order_dt) * 0.16))
figure5_write_tsv(dt, file.path(paths$metadata, "figure5_six_panel_C_heatmap_matrix.tsv.gz"))
figure5_write_tsv(order_dt[, .(entity, display_order, axis, entity_type)],
                  file.path(paths$metadata, "figure5_six_panel_C_programme_order.tsv"))
figure5_write_json(list(
  panel = "5C",
  title = "Temporal organization of frozen TF, regulon and target-gene programmes",
  frozen_signatures = TRUE,
  source_namespace = "figure5_temporal_positioning_onset_fix",
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_C_report.json"))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5_panel_C.rds"))

