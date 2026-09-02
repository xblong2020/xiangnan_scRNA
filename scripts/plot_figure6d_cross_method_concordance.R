#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

x <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6d_cross_method_concordance.tsv"))
p <- ggplot(x, aes(top50_jaccard, spearman_rank_correlation, colour = axis, size = top50_overlap, label = tf)) +
  geom_hline(yintercept = 0, linewidth = .35, linetype = 2, colour = neutral_gray) +
  geom_point(alpha = .85) +
  ggrepel::geom_text_repel(size = 2.45, colour = dark_text, max.overlaps = Inf, box.padding = .25, seed = 20260805) +
  scale_colour_manual(values = c(axis_palette, control = neutral_gray), name = "Perturbation class") +
  scale_size_continuous(range = c(1.8, 5), name = "Top-50 overlap") +
  labs(title = "Cross-model perturbation concordance", x = "Top-50 gene-set Jaccard",
    y = "Rank correlation of perturbation magnitude",
    caption = "CellOracle is signed; scTenifoldKnk manifold distance is unsigned. Sign concordance is therefore not estimable.") +
  figure6_theme() + theme(legend.position = "right")
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6d_cross_method_concordance")
figure6_save(p, out_dir, "figure6d_cross_method_concordance", 6.4, 4.7)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6d_plot.rds"))

