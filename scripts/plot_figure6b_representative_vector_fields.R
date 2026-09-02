#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

cells <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_cells.tsv.gz"))
grid <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_grid.tsv.gz"))
cells[, celloracle_state := factor(celloracle_state, levels = names(state_palette))]
titles <- c(identity_axis = "HNF4A knockout", stress_axis = "EGR1 knockout", sox4_axis = "SOX4 knockout")
cells[, panel := factor(titles[axis], levels = unname(titles))]
grid[, panel := factor(titles[axis], levels = unname(titles))]
xlim <- range(cells$umap_1, finite = TRUE); ylim <- range(cells$umap_2, finite = TRUE)
p <- ggplot() +
  geom_point(data = cells, aes(umap_1, umap_2, colour = celloracle_state), size = 0.28, alpha = 0.35, stroke = 0) +
  geom_segment(data = grid[show %in% TRUE], aes(x = grid_x, y = grid_y, xend = plot_xend, yend = plot_yend),
    arrow = grid::arrow(length = grid::unit(0.037, "inches"), type = "open"),
    linewidth = 0.32, colour = dark_text, lineend = "round") +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = state_palette, drop = FALSE, name = "Baseline state") +
  coord_fixed(xlim = xlim, ylim = ylim, clip = "off") +
  labs(title = "Representative CellOracle knockout vector fields", x = "CellOracle UMAP 1", y = "CellOracle UMAP 2",
    caption = "Frozen, previously audited Figure 2–4 fields; common coordinate range and arrow-length scaling.") +
  figure6_theme() +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(size = 8.3, colour = dark_text)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.2, alpha = 1)))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6b_representative_vector_fields")
figure6_save(p, out_dir, "figure6b_representative_vector_fields", 8.4, 3.3)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6b_plot.rds"))
