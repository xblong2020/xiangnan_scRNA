#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(ggplot2); library(ggsci); library(scales)})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))
root <- figure2_project_root()
input_path <- figure2_get_arg("--input", file.path(root,
  "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_perturbation_genes.tsv"))
sensitivity_path <- figure2_get_arg("--sensitivity-input", file.path(root,
  "metadata/driver/sctenifoldknk_module7_2_malignant_like_perturbation_genes.tsv"))
out_dir <- figure2_get_arg("--out-dir", file.path(root, "metadata/driver/figure2e_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2e_hnf4a"))
sens_out <- figure2_get_arg("--sensitivity-out-dir",
  file.path(root, "metadata/driver/figure2e_hnf4a_sensitivity"))
sens_fig <- figure2_get_arg("--sensitivity-figure-dir",
  file.path(root, "figures/driver/figure2e_hnf4a_sensitivity"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")
top_n <- as.integer(figure2_get_arg("--top-n", "20"))
fdr <- as.numeric(figure2_get_arg("--fdr-cutoff", "0.05"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_out, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_fig, recursive = TRUE, showWarnings = FALSE)

normalize_columns <- function(dat) {
  if ("Distance" %in% names(dat) && !"distance" %in% names(dat)) names(dat)[names(dat) == "Distance"] <- "distance"
  if ("p_adj" %in% names(dat) && !"p.adj" %in% names(dat)) names(dat)[names(dat) == "p_adj"] <- "p.adj"
  dat$distance <- as.numeric(dat$distance); dat$p.adj <- as.numeric(dat$p.adj)
  dat
}
select_sig <- function(path) {
  dat <- normalize_columns(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  required <- c("tf", "gene", "distance", "p.adj", "subset")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x <- dat[dat$tf == target_tf & dat$gene != target_tf &
             is.finite(dat$distance) & is.finite(dat$p.adj) & dat$p.adj < fdr, , drop = FALSE]
  x[order(x$distance, decreasing = TRUE, x$gene), , drop = FALSE]
}
sig <- select_sig(input_path)
sens <- select_sig(sensitivity_path)
plot_dat <- head(sig, top_n)
if (nrow(plot_dat)) {
  plot_dat$rank <- seq_len(nrow(plot_dat))
  plot_dat$minus_log10_fdr <- -log10(pmax(plot_dat$p.adj, .Machine$double.xmin))
  plot_dat$significance <- ifelse(plot_dat$p.adj < .001, "***",
    ifelse(plot_dat$p.adj < .01, "**", "*"))
  plot_dat$gene_plot <- factor(plot_dat$gene, levels = rev(plot_dat$gene))
}
source_path <- file.path(out_dir, "figure2e_hnf4a_significant_perturbed_genes.tsv")
write.table(plot_dat, source_path, sep = "\t", quote = FALSE, row.names = FALSE)

stem <- file.path(figure_dir, "figure2e_hnf4a_significant_perturbed_genes")
figure_generated <- FALSE
title <- if (nrow(sig) >= top_n) "Top 20 significantly perturbed genes" else
  paste0("Significantly perturbed genes after HNF4A knockout (n = ", nrow(sig), ")")
if (nrow(plot_dat)) {
  lancet <- ggsci::pal_lancet("lanonc")(9)
  cols <- c(lancet[1], lancet[4], lancet[3], lancet[6], lancet[2])
  x_max <- max(plot_dat$distance)
  height <- max(2.6, 1.3 + 0.195 * nrow(plot_dat))
  p <- ggplot(plot_dat, aes(distance, gene_plot)) +
    geom_col(aes(fill = distance), width = .68, alpha = .90) +
    geom_point(aes(fill = distance), shape = 21, size = 2.25, stroke = .45, colour = "white") +
    geom_text(aes(label = significance), x = plot_dat$distance + x_max * .025,
              hjust = 0, size = 2.8, colour = "grey20") +
    scale_x_continuous(labels = label_scientific(digits = 2),
                       expand = expansion(mult = c(0, .16))) +
    scale_fill_gradientn(colours = cols, name = "Distance") +
    labs(x = "Manifold alignment distance", y = NULL, title = title, tag = "Figure 2E") +
    coord_cartesian(clip = "off") + figure2_theme() +
    theme(axis.text.y = element_text(size = 7.7), plot.margin = margin(6, 15, 6, 6))
  figure2_save(p, figure_dir, "figure2e_hnf4a_significant_perturbed_genes", 5.4, height, tiff = TRUE)
  figure_generated <- TRUE
}

main_genes <- unique(sig$gene); sens_genes <- unique(sens$gene)
union_n <- length(union(main_genes, sens_genes))
jaccard <- if (union_n) length(intersect(main_genes, sens_genes)) / union_n else NA_real_
overlap <- intersect(main_genes, sens_genes)
rank_cor <- NA_real_
if (length(overlap) >= 3) {
  rank_cor <- cor(match(overlap, main_genes), match(overlap, sens_genes), method = "spearman")
}
comparison <- data.frame(
  target_tf = target_tf, main_subset = "normal_reference",
  sensitivity_subset = "malignant_like", n_significant_main = length(main_genes),
  n_significant_sensitivity = length(sens_genes), n_overlap = length(overlap),
  jaccard_overlap = jaccard, distance_rank_spearman = rank_cor
)
write.table(comparison, file.path(sens_out, "figure2e_hnf4a_sensitivity_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
counts <- data.frame(subset = c("normal_reference", "malignant_like"),
                     n_significant = c(length(main_genes), length(sens_genes)))
ps <- ggplot(counts, aes(subset, n_significant, fill = subset)) +
  geom_col(width = .62) + geom_text(aes(label = n_significant), vjust = -.25, size = 3) +
  scale_fill_manual(values = c(normal_reference = "#56B4E9", malignant_like = "#D55E00"),
                    guide = "none") +
  labs(x = NULL, y = "FDR-significant genes", title = "HNF4A network sensitivity") +
  figure2_theme()
figure2_save(ps, sens_fig, "figure2e_hnf4a_sensitivity_significant_gene_counts", 4.2, 3.4)

report <- list(
  module = "Figure 2E", target_tf = target_tf,
  method = "scTenifoldKnk HNF4A virtual knockout network perturbation evidence",
  input = figure2_norm_path(input_path), subset = "normal_reference",
  target_gene_excluded = TRUE, fdr_cutoff = fdr,
  selection = "Only p.adj < 0.05 genes, descending manifold alignment distance",
  n_significant_excluding_target = nrow(sig), n_plotted = nrow(plot_dat),
  figure_generated = figure_generated, dynamic_title = title,
  top_genes = as.list(plot_dat$gene),
  sensitivity = as.list(comparison[1, ]),
  source_table = figure2_norm_path(source_path),
  caveat = "Manifold distance is non-directional; results support network displacement, not activation or suppression."
)
figure2_write_json(report, file.path(out_dir, "figure2e_hnf4a_report.json"))
figure2_write_json(list(module = "Figure 2E sensitivity", target_tf = target_tf,
                        comparison = as.list(comparison[1, ]),
                        caveat = "State-specific scTenifoldKnk virtual knockout comparison."),
                   file.path(sens_out, "figure2e_hnf4a_sensitivity_report.json"))
