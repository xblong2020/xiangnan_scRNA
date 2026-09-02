#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages({library(ComplexUpset); library(patchwork)})

sets <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6e_gene_sets.tsv"))
jac <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6e_jaccard_matrix.tsv"))
path_jac <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6e_pathway_similarity.tsv"))
focus <- c("HNF4A", "PPARA", "EGR1", "CEBPB", "AP-1 aggregate", "SOX4")
wide <- dcast(unique(sets[, .(gene, perturbation)])[, value := TRUE], gene ~ perturbation, value.var = "value", fill = FALSE)
for (z in focus) if (!z %in% names(wide)) wide[, (z) := FALSE]
p_upset <- ComplexUpset::upset(
  as.data.frame(wide), focus, min_size = 1, width_ratio = .25, n_intersections = 18,
  sort_intersections_by = "cardinality",
  themes = ComplexUpset::upset_modify_themes(list(
    intersections_matrix = theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  )),
  base_annotations = list("Intersection size" = ComplexUpset::intersection_size(text = list(size = 2.5))),
  set_sizes = ComplexUpset::upset_set_size() + scale_y_continuous(expand = expansion(mult = c(.02, .12)))
) + labs(title = "Unified-background perturbation gene sets") + figure6_theme()
jac[, `:=`(set_1 = factor(set_1, levels = focus), set_2 = factor(set_2, levels = rev(focus)))]
p_heat <- ggplot(jac, aes(set_1, set_2, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = .45) + geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.35, colour = dark_text) +
  scale_fill_gradient(low = "white", high = auxiliary_purple, limits = c(0, 1), name = "Jaccard") +
  labs(title = "Gene-set similarity", x = NULL, y = NULL,
    caption = if (any(path_jac$status == "estimable")) "FDR pathway-set similarity is available in the panel source data." else "Pathway similarity: not estimable (no FDR-significant pathways).") +
  figure6_theme() + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
p <- wrap_elements(full = p_upset) + p_heat + plot_layout(widths = c(1.45, 1)) +
  plot_annotation(title = "Perturbation target and pathway overlap", theme = theme(plot.title = element_text(size = 10, hjust = .5)))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6e_target_pathway_overlap")
figure6_save(p, out_dir, "figure6e_target_pathway_overlap", 9.1, 4.9)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6e_plot.rds"))
