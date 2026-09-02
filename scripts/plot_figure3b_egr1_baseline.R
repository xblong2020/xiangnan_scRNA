#!/usr/bin/env Rscript

## Figure 3B: exact common baseline developmental field reused for EGR1.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "scripts", "figure3_egr1_common.R"))

data_dir <- normalizePath(
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3b_egr1")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3b_egr1")),
  mustWork = FALSE
)
seed <- as.integer(figure3_get_arg("--seed", "15071990"))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cells_path <- file.path(data_dir, "figure3b_egr1_plot_cells.tsv.gz")
umap_grid_path <- file.path(data_dir, "figure3b_egr1_baseline_grid_umap.tsv.gz")
tsne_grid_path <- file.path(data_dir, "figure3b_egr1_baseline_grid_tsne.tsv.gz")
if (!all(file.exists(c(cells_path, umap_grid_path, tsne_grid_path)))) {
  stop("Prepared Figure 3B common-baseline tables are missing.")
}
cells <- read.delim(gzfile(cells_path), stringsAsFactors = FALSE, check.names = FALSE)
umap_grid <- read.delim(gzfile(umap_grid_path), stringsAsFactors = FALSE)
tsne_grid <- read.delim(gzfile(tsne_grid_path), stringsAsFactors = FALSE)
cells$state_label <- factor(
  unname(figure3_state_labels[cells$celloracle_state]),
  levels = unname(figure3_state_labels[figure3_state_order])
)
palette_labelled <- setNames(
  unname(figure3_state_palette),
  unname(figure3_state_labels[figure3_state_order])
)

plot_field <- function(df, x_col, y_col, grid, x_lab, tag = NULL) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$keep %in% TRUE, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = state_label)) +
    geom_point(size = 0.48, alpha = 0.48, stroke = 0) +
    geom_segment(
      data = arrows,
      aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE,
      colour = "grey15",
      linewidth = 0.28,
      arrow = arrow(length = grid::unit(0.05, "inches"), type = "closed")
    ) +
    ggrepel::geom_text_repel(
      data = centres,
      aes(x = label_x, y = label_y, label = state),
      inherit.aes = FALSE,
      size = 2.4,
      family = "sans",
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = Inf,
      seed = seed,
      colour = "grey15"
    ) +
    scale_colour_manual(values = palette_labelled, drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = "Baseline developmental field", tag = tag) +
    figure3_theme() +
    theme(legend.position = "right")
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "Figure 3B")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1")
figure3_save(p_umap, figure_dir, "figure3b_egr1_baseline_umap", 5.7, 4.5)
figure3_save(p_tsne, figure_dir, "figure3b_egr1_baseline_tsne", 5.7, 4.5)
p_combined <- p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right")
figure3_save(p_combined, figure_dir, "figure3b_egr1_baseline_umap_tsne", 10.8, 4.5)

report <- list(
  module = "Figure 3B",
  target_tf = "EGR1",
  title = "Baseline developmental field",
  plotting_language = "R",
  r_version = R.version.string,
  n_cells = nrow(cells),
  source_contract = "Exact value-level reuse of validated common baseline cells, UMAP, t-SNE, pseudotime, grids, arrows, and masks",
  vector_field = list(
    n_grid = 20,
    k_neighbors = 50,
    density_quantile = 0.70,
    seed = seed,
    baseline_tf_independent = TRUE
  ),
  state_palette = as.list(figure3_state_palette),
  outputs = list(
    umap_pdf = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap.pdf")),
    umap_png = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap.png")),
    umap_svg = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap.svg")),
    tsne_pdf = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_tsne.pdf")),
    tsne_png = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_tsne.png")),
    tsne_svg = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_tsne.svg")),
    combined_pdf = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap_tsne.pdf")),
    combined_png = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap_tsne.png")),
    combined_svg = figure3_norm_path(file.path(figure_dir, "figure3b_egr1_baseline_umap_tsne.svg"))
  ),
  caveat = "Baseline developmental field is TF-independent; the t-SNE coordinates were reused and not recalculated."
)
figure3_write_json(report, file.path(data_dir, "figure3b_egr1_r_plot_report.json"))
message("Figure 3B written to: ", figure_dir)

