#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
paths <- figure5_onset_fix_paths(root)
dt <- fread(file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))
dt[, landmark_label := factor(landmark, levels = c("onset", "maximum_slope", "t50"), labels = c("Onset", "Maximum slope (primary)", "t50"))]
dt[, comparison := factor(comparison, levels = c("A before B", "B before C", "A before C"))]
dt[, label := ifelse(is.finite(probability), sprintf("P=%.2f\nΔ[%.2f, %.2f]\ntie=%.1f%%", probability, delta_q025, delta_q975, 100 * tie_fraction), "NA")]

p <- ggplot(dt, aes(comparison, landmark_label, fill = evidence_grade)) +
  geom_tile(colour = lancet_palette[8], linewidth = 0.45) + geom_text(aes(label = label), size = 2.45, lineheight = 0.9) +
  scale_fill_manual(values = evidence_palette, drop = FALSE) +
  labs(title = "E  Tie-aware temporal-precedence matrix", subtitle = "Cell values: precedence probability, bootstrap delta quantiles, and tie fraction",
       x = NULL, y = NULL, fill = "Evidence") + theme_figure5() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

out_dir <- file.path(paths$figures, "figure5e_precedence_matrix")
export_figure5_plot(p, file.path(out_dir, "figure5e_precedence_matrix"), 5.5, 3.5)
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5e.rds"))
