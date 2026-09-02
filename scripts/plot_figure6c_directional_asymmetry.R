#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages(library(patchwork))

x <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry.tsv"))
x[, comparison := factor(comparison, levels = rev(comparison))]
pair_long <- rbindlist(list(
  x[, .(comparison, direction = "Forward", absolute_effect = forward_absolute_effect, source_axis)],
  x[, .(comparison, direction = "Reverse", absolute_effect = reverse_absolute_effect, source_axis)]
))
p1 <- ggplot(x, aes(y = comparison)) +
  geom_segment(aes(x = reverse_absolute_effect, xend = forward_absolute_effect, yend = comparison), colour = neutral_gray, linewidth = 0.7) +
  geom_point(aes(x = forward_absolute_effect, colour = source_axis), size = 2.8) +
  geom_point(aes(x = reverse_absolute_effect), colour = neutral_gray, size = 2.8, shape = 21, fill = "white") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  labs(title = "Forward versus reverse magnitude", x = "Absolute standardized effect", y = NULL,
    caption = "Filled: forward; open: reverse") + figure6_theme()
p2 <- ggplot(x, aes(directional_asymmetry_score, comparison, colour = classification)) +
  geom_vline(xintercept = 0, linetype = 2, colour = neutral_gray, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.13, linewidth = 0.6) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = c("Forward-dominant" = unname(axis_palette["stress_axis"]), "Reverse-dominant" = unname(axis_palette["identity_axis"]),
    "Symmetric/unresolved" = unname(neutral_gray)), drop = FALSE) +
  labs(title = "Directional asymmetry score", x = "|Forward| − |Reverse| (95% CI)", y = NULL, colour = NULL,
    caption = "Sample bootstrap stratified by dataset") + figure6_theme() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
p <- p1 + p2 + plot_layout(widths = c(1.05, 1)) + plot_annotation(title = "Forward–reverse directional asymmetry", theme = theme(plot.title = element_text(size = 10, hjust = .5)))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6c_directional_asymmetry")
figure6_save(p, out_dir, "figure6c_directional_asymmetry", 8.2, 3.5)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6c_plot.rds"))
