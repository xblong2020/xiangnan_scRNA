#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
paths <- figure5_onset_fix_paths(root)
dt <- fread(file.path(paths$metadata, "figure5f_method_concordance.tsv"))
method_order <- unique(dt$method)
dt[, method := factor(method, levels = rev(method_order))]
dt[, comparison := factor(comparison, levels = c("A before B", "B before C", "A before C"))]
dt[, label := ifelse(is.finite(probability), sprintf("%.2f\nn=%d", probability, n_valid), "NA")]

p <- ggplot(dt, aes(comparison, method, fill = status)) +
  geom_tile(colour = lancet_palette[8], linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.15, lineheight = 0.9) +
  scale_fill_manual(values = evidence_palette, drop = FALSE) +
  labs(
    title = "F  Corrected method and sensitivity concordance",
    subtitle = "All fitted GAMs use REML; n is the valid fitted-unit/resample count and small-n rows are descriptive",
    x = NULL, y = NULL, fill = "Evidence"
  ) +
  theme_figure5(7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

out_dir <- file.path(paths$figures, "figure5f_method_concordance")
export_figure5_plot(p, file.path(out_dir, "figure5f_method_concordance"), 6.4, max(4.8, uniqueN(dt$method) * 0.25))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5f.rds"))
