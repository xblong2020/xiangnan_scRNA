#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(patchwork))
paths <- figure5_paths(root)

states <- data.table(
  state = c("Normal/\nreference", "Stressed/\ninjured", "Regenerative/\nprogenitor", "Proliferating\ncandidate", "Malignant/\nmalignant-like"),
  xmin = seq(0, 0.8, by = 0.2), xmax = seq(0.2, 1, by = 0.2), y = 3.55
)
bands <- data.table(
  axis = factor(c("identity_loss", "stress_transition", "sox4_stabilization"), levels = names(axis_palette)),
  label = c("Axis A  HNF4A/PPARA identity loss", "Axis B  AP-1/CEBPB/EGR1 stress transition", "Axis C  SOX4 malignant-state stabilization"),
  xmin = c(0.08, 0.20, 0.50), xmax = c(0.72, 0.86, 0.98), y = c(2.65, 2.02, 1.39)
)
evidence <- c("Consensus pseudotime", "DPT (not available)", "Monocle3", "Slingshot", "CellRank fate",
              "CytoTRACE2", "CNV-supported fate", "Patient pseudobulk", "Bootstrap ordering", "Leave-one-dataset-out")

p <- ggplot() +
  geom_rect(data = states, aes(xmin = xmin, xmax = xmax, ymin = y - 0.24, ymax = y + 0.24),
            fill = alpha(lancet_palette[8], 0.18), colour = lancet_palette[8], linewidth = 0.35) +
  geom_text(data = states, aes(x = (xmin + xmax) / 2, y = y, label = state), size = 2.55, lineheight = 0.95) +
  geom_segment(data = states[-nrow(states)], aes(x = xmax, xend = xmax + 0.035, y = y, yend = y),
               arrow = grid::arrow(length = grid::unit(1.5, "mm")), linewidth = 0.35, colour = lancet_palette[9]) +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = y - 0.22, ymax = y + 0.22, fill = axis),
            alpha = 0.38, colour = NA) +
  geom_segment(data = bands, aes(x = xmin, xend = xmax, y = y, yend = y, colour = axis), linewidth = 1.1, alpha = 0.75) +
  geom_text(data = bands, aes(x = xmin + 0.015, y = y, label = label), hjust = 0, size = 2.65, fontface = "bold") +
  annotate("rect", xmin = 1.04, xmax = 1.40, ymin = 0.75, ymax = 3.90, fill = alpha(lancet_palette[6], 0.12), colour = lancet_palette[8], linewidth = 0.4) +
  annotate("text", x = 1.22, y = 3.68, label = "Triangulated evidence", fontface = "bold", size = 2.8) +
  annotate("text", x = 1.07, y = seq(3.38, 1.10, length.out = length(evidence)), label = paste0("\u2022 ", evidence), hjust = 0, size = 2.35) +
  annotate("text", x = 0.70, y = 0.35, label = "Pseudotemporal positioning does not establish physical time or direct causality.",
           size = 2.55, fontface = "italic", colour = lancet_palette[9]) +
  scale_fill_manual(values = axis_palette, guide = "none") +
  scale_colour_manual(values = axis_palette, guide = "none") +
  coord_cartesian(xlim = c(0, 1.42), ylim = c(0.15, 4.05), clip = "off") +
  labs(title = "A  Temporal-positioning analysis framework", x = NULL, y = NULL) +
  theme_void(base_family = "sans", base_size = 7.5) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0), plot.margin = margin(5, 8, 5, 8))

out_dir <- file.path(paths$figures, "figure5a_temporal_framework")
outputs <- export_figure5_plot(p, file.path(out_dir, "figure5a_temporal_framework"), 7.2, 4.1)
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5a.rds"))
figure5_write_json(list(panel = "5A", schematic = TRUE, axis_ranges_overlap = TRUE, outputs = as.list(outputs)),
                   file.path(paths$metadata, "figure5a_temporal_framework_report.json"))
