#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))
root <- figure2_project_root()
data_dir <- figure2_get_arg("--data-dir", file.path(root, "metadata/driver/figure2c_hnf4a"))
b_dir <- figure2_get_arg("--figure2b-data-dir", file.path(root, "metadata/driver/figure2b_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2c_hnf4a"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")

cells <- read.delim(gzfile(file.path(data_dir, "figure2c_hnf4a_cells_with_tsne_projection.tsv.gz")),
                    stringsAsFactors = FALSE)
ug <- read.delim(gzfile(file.path(data_dir, "figure2c_hnf4a_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tg <- read.delim(gzfile(file.path(data_dir, "figure2c_hnf4a_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)
ur <- read.delim(gzfile(file.path(b_dir, "figure2b_hnf4a_baseline_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tr <- read.delim(gzfile(file.path(b_dir, "figure2b_hnf4a_baseline_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)

score_grid <- function(perturb, reference) {
  if (nrow(perturb) != nrow(reference) ||
      any(perturb$grid_x != reference$grid_x) || any(perturb$grid_y != reference$grid_y)) {
    stop("Perturbation and baseline grids are not pointwise aligned")
  }
  out <- perturb
  out$inner_product_score_grid <- perturb$flow_x * reference$unit_x +
    perturb$flow_y * reference$unit_y
  out$keep_score <- (perturb$keep %in% TRUE) & is.finite(out$inner_product_score_grid)
  out$show_score <- out$keep_score & (perturb$show %in% TRUE)
  out
}
us <- score_grid(ug, ur); ts <- score_grid(tg, tr)
all_scores <- c(us$inner_product_score_grid[us$keep_score], ts$inner_product_score_grid[ts$keep_score])
score_limit <- max(abs(all_scores), na.rm = TRUE)
if (!is.finite(score_limit) || score_limit == 0) score_limit <- 1
write.table(us, gzfile(file.path(data_dir, "figure2c_hnf4a_inner_product_grid_umap.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ts, gzfile(file.path(data_dir, "figure2c_hnf4a_inner_product_grid_tsne.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_score <- function(df, x_col, y_col, grid, x_lab, tag = NULL) {
  points <- grid[grid$show_score %in% TRUE, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(colour = "grey82", size = 0.45, alpha = 0.22) +
    geom_point(data = points, aes(x = grid_x, y = grid_y, colour = inner_product_score_grid),
               inherit.aes = FALSE, size = 2, alpha = 0.95) +
    scale_colour_gradient2(low = figure2_diverging["low"], mid = figure2_diverging["mid"],
      high = figure2_diverging["high"], midpoint = 0, limits = c(-score_limit, score_limit),
      oob = scales::squish, name = "Inner product score (PS)") +
    geom_segment(data = points, aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE, colour = "black", linewidth = 0.28, alpha = 0.75,
      arrow = arrow(length = grid::unit(0.045, "inches"), type = "closed")) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = "HNF4A perturbation score", tag = tag) +
    figure2_theme() + theme(legend.position = "right")
}

p_umap <- plot_score(cells, "umap_1", "umap_2", us, "CellOracle UMAP1", "Figure 2C")
p_tsne <- plot_score(cells, "tsne_1", "tsne_2", ts, "Expression t-SNE1")
figure2_save(p_umap, figure_dir, "figure2c_hnf4a_inner_product_umap", 5.7, 4.5)
figure2_save(p_tsne, figure_dir, "figure2c_hnf4a_inner_product_tsne", 5.7, 4.5)
figure2_save(p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right"),
             figure_dir, "figure2c_hnf4a_inner_product_umap_tsne", 10.8, 4.5)

uv <- us$inner_product_score_grid[us$keep_score]
tv <- ts$inner_product_score_grid[ts$keep_score]
report <- list(
  module = "Figure 2C score-coloured companion", target_tf = target_tf,
  definition = "HNF4A perturbation vector dot baseline developmental vector",
  n_cells = nrow(cells), n_grid_points = nrow(us),
  score_range = range(uv), score_range_umap = range(uv), score_range_tsne = range(tv),
  positive_ps_fraction = mean(uv > 0), negative_ps_fraction = mean(uv < 0),
  zero_center = 0, colour_scale = as.list(figure2_diverging),
  umap_method = "Saved HNF4A CellOracle delta_embedding aggregated to matched baseline grid",
  tsne_method = "Supplementary local UMAP-to-t-SNE Jacobian projection; not native simulation",
  caveat = "PS is a computationally inferred directional alignment score from virtual knockout."
)
figure2_write_json(report, file.path(data_dir, "figure2c_hnf4a_inner_product_report.json"))
