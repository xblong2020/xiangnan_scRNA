#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))

root <- figure2_project_root()
data_dir <- figure2_get_arg("--data-dir", file.path(root, "metadata/driver/figure2b_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2b_hnf4a"))
seed <- as.integer(figure2_get_arg("--seed", "15071990"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim(gzfile(file.path(data_dir, "figure2b_hnf4a_plot_cells.tsv.gz")),
                    stringsAsFactors = FALSE, check.names = FALSE)
umap_grid <- read.delim(gzfile(file.path(data_dir, "figure2b_hnf4a_baseline_grid_umap.tsv.gz")),
                        stringsAsFactors = FALSE)
tsne_grid <- read.delim(gzfile(file.path(data_dir, "figure2b_hnf4a_baseline_grid_tsne.tsv.gz")),
                        stringsAsFactors = FALSE)
cells$state_label <- factor(unname(figure2_state_labels[cells$celloracle_state]),
                            levels = unname(figure2_state_labels[figure2_state_order]))
palette_labelled <- setNames(unname(figure2_state_palette),
                             unname(figure2_state_labels[figure2_state_order]))

plot_field <- function(df, x_col, y_col, grid, x_lab, tag = NULL) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$keep %in% TRUE, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = state_label)) +
    geom_point(size = 0.48, alpha = 0.48, stroke = 0) +
    geom_segment(data = arrows,
      aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE, colour = "grey15", linewidth = 0.28,
      arrow = arrow(length = grid::unit(0.05, "inches"), type = "closed")) +
    ggrepel::geom_text_repel(data = centres,
      aes(x = label_x, y = label_y, label = state), inherit.aes = FALSE,
      size = 2.4, family = "sans", box.padding = 0.25, point.padding = 0.15,
      min.segment.length = Inf, seed = seed, colour = "grey15") +
    scale_colour_manual(values = palette_labelled, drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = "Baseline developmental field", tag = tag) +
    figure2_theme() + theme(legend.position = "right")
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "Figure 2B")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1")
figure2_save(p_umap, figure_dir, "figure2b_hnf4a_baseline_umap", 5.7, 4.5)
figure2_save(p_tsne, figure_dir, "figure2b_hnf4a_baseline_tsne", 5.7, 4.5)
p_combined <- p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right")
figure2_save(p_combined, figure_dir, "figure2b_hnf4a_baseline_umap_tsne", 10.8, 4.5)

report <- list(
  module = "Figure 2B", target_tf = target_tf, title = "Baseline developmental field",
  plotting_language = "R", r_version = R.version.string, n_cells = nrow(cells),
  source_contract = "Exact value-level reuse of SOX4 Figure 2B cells, coordinates, pseudotime and grids",
  vector_field = list(n_grid = 20, k_neighbors = 50, density_quantile = 0.70,
                      seed = seed, baseline_tf_independent = TRUE),
  state_palette = as.list(figure2_state_palette),
  outputs = list(
    umap_pdf = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap.pdf")),
    umap_png = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap.png")),
    umap_svg = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap.svg")),
    tsne_pdf = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_tsne.pdf")),
    tsne_png = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_tsne.png")),
    tsne_svg = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_tsne.svg")),
    combined_pdf = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap_tsne.pdf")),
    combined_png = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap_tsne.png")),
    combined_svg = figure2_norm_path(file.path(figure_dir, "figure2b_hnf4a_baseline_umap_tsne.svg"))
  ),
  caveat = "Baseline developmental field is TF-independent."
)
figure2_write_json(report, file.path(data_dir, "figure2b_hnf4a_r_plot_report.json"))
