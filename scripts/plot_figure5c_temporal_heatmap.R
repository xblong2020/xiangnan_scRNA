#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
paths <- figure5_onset_fix_paths(root)
dt <- as.data.table(readRDS(file.path(paths$processed, "figure5c_heatmap_predictions.rds")))
order_dt <- unique(dt[, .(entity, display_order, axis, entity_type)])[order(display_order)]
dt[, entity_factor := factor(entity, levels = rev(order_dt$entity))]
order_dt[, entity_factor := factor(entity, levels = rev(order_dt$entity))]

p <- ggplot(dt, aes(pseudotime, entity_factor, fill = row_z)) +
  geom_tile() +
  geom_segment(data = order_dt, aes(x = -0.035, xend = -0.012, y = entity_factor,
                                    yend = entity_factor, colour = axis),
               inherit.aes = FALSE, linewidth = 2.1) +
  scale_fill_gradientn(colours = figure5_diverging_palette, limits = c(-2.5, 2.5), oob = squish, name = "Row z-score") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  coord_cartesian(xlim = c(-0.04, 1), clip = "off") +
  labs(title = "C  Temporal ordering of frozen TFs and targets", x = "Normalized pseudotime", y = NULL) +
  theme_figure5(7) + theme(axis.text.y = element_text(size = 5.8), plot.margin = margin(5, 6, 5, 18))

out_dir <- file.path(paths$figures, "figure5c_temporal_heatmap")
export_figure5_plot(p, file.path(out_dir, "figure5c_temporal_heatmap"), 5.3, max(5.8, nrow(order_dt) * 0.16))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5c.rds"))
