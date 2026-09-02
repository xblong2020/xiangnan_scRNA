#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

x <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6h_comparison_table.tsv"))
x[, y := .N:1]
x[, published_label := vapply(published_foxm1_cebpb, function(z) paste(strwrap(z, width = 34), collapse = "\n"), character(1))]
x[, current_label := vapply(current_three_axis, function(z) paste(strwrap(z, width = 38), collapse = "\n"), character(1))]
rel_cols <- c("complementary" = auxiliary_cyan, "partially overlapping" = auxiliary_purple,
  "distinct analytical scope" = axis_palette["sox4_axis"], "distinct evidence level" = neutral_gray, "unresolved" = lancet_palette[6])
p <- ggplot(x, aes(y = y)) +
  geom_label(aes(x = .9, label = published_label), fill = auxiliary_purple, colour = "white", size = 2.1,
    label.size = .2, label.padding = grid::unit(.12, "lines"), lineheight = .9) +
  geom_segment(aes(x = 1.18, xend = 1.92, yend = y, colour = relationship), linewidth = .8) +
  geom_point(aes(x = 1.55, colour = relationship), size = 2.3) +
  geom_label(aes(x = 2.2, label = current_label), fill = axis_palette["stress_axis"], colour = "white", size = 2.1,
    label.size = .2, label.padding = grid::unit(.12, "lines"), lineheight = .9) +
  geom_text(aes(x = .32, label = dimension), hjust = 1, size = 2.35, colour = dark_text) +
  scale_colour_manual(values = rel_cols, name = "Relationship") +
  scale_x_continuous(limits = c(.25,2.85), breaks = c(.9,2.2), labels = c("Published FOXM1/CEBPB axis", "Current three-axis model"), position = "top") +
  scale_y_continuous(expand = expansion(mult = c(.02,.05))) +
  coord_cartesian(clip = "off") +
  labs(title = "Comparison with the published FOXM1/CEBPB plasticity axis",
    caption = "Published model: Zhang et al., J Hepatol 2026, PMID 41043722. Current architecture is computationally inferred; FOXM1 was not perturbed here.",
    x = NULL, y = NULL) +
  theme_void(base_family = "sans") + theme(plot.title=element_text(size=10,hjust=.5,colour=dark_text), plot.caption=element_text(size=6.8,hjust=0,colour=dark_text),
    axis.text.x.top=element_text(size=8,colour=dark_text), legend.position="bottom", legend.text=element_text(size=7), legend.title=element_text(size=7.2), plot.margin=margin(8,14,8,86))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6h_foxm1_cebpb_comparison")
figure6_save(p, out_dir, "figure6h_foxm1_cebpb_comparison", 9.8, 7.4)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6h_plot.rds"))
