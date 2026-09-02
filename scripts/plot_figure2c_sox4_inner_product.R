#!/usr/bin/env Rscript

## Score-coloured Figure 2C companion.
## The continuous colour represents the CellOracle perturbation score (PS),
## i.e. the inner product between the perturbation and developmental vectors.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(), mustWork = FALSE)
if (!dir.exists(file.path(PROJECT_ROOT, "scripts"))) PROJECT_ROOT <- normalizePath(getwd(), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

data_dir <- normalizePath(get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2c_sox4")), mustWork = FALSE)
b_data_dir <- normalizePath(get_arg("--figure2b-data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2b_sox4")), mustWork = FALSE)
figure_dir <- normalizePath(get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure2c_sox4")), mustWork = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
cells <- read.delim(gzfile(file.path(data_dir, "figure2c_sox4_matched_cells.tsv.gz")), stringsAsFactors = FALSE)
umap_grid <- read.delim(gzfile(file.path(data_dir, "figure2c_sox4_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tsne_grid <- read.delim(gzfile(file.path(data_dir, "figure2c_sox4_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)
umap_ref <- read.delim(gzfile(file.path(b_data_dir, "figure2b_sox4_baseline_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tsne_ref <- read.delim(gzfile(file.path(b_data_dir, "figure2b_sox4_baseline_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)

## Project's existing state colours are retained only for the faint background;
## the foreground colour encodes PS on a continuous coolwarm-like scale.
state_palette <- c(
  normal_reference = "#B8B8B8",
  stressed_injured = "#56B4E9",
  regenerative_progenitor = "#009E73",
  proliferating_candidate = "#E69F00",
  malignant_or_malignant_like = "#D55E00"
)

## The native CellOracle definition is dot(perturbation flow, reference flow).
## Our saved R baseline field stores a normalized reference direction, so this
## is a direction-and-magnitude PS approximation on the matched grid.
score_grid <- function(perturb, reference) {
  out <- perturb
  out$inner_product_score_grid <- perturb$flow_x * reference$unit_x + perturb$flow_y * reference$unit_y
  out$keep_score <- (perturb$keep %in% TRUE) & is.finite(out$inner_product_score_grid)
  out$show_score <- out$keep_score & (perturb$show %in% TRUE)
  out
}
umap_score <- score_grid(umap_grid, umap_ref)
tsne_score <- score_grid(tsne_grid, tsne_ref)
all_scores <- c(umap_score$inner_product_score_grid[umap_score$keep_score], tsne_score$inner_product_score_grid[tsne_score$keep_score])
score_limit <- max(abs(all_scores), na.rm = TRUE)
if (!is.finite(score_limit) || score_limit == 0) score_limit <- 1
write.table(umap_score, gzfile(file.path(data_dir, "figure2c_sox4_inner_product_grid_umap.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne_score, gzfile(file.path(data_dir, "figure2c_sox4_inner_product_grid_tsne.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

plot_score <- function(df, x_col, y_col, grid, x_label, title) {
  points <- grid[grid$show_score, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(colour = "grey82", size = 0.45, alpha = 0.22, show.legend = FALSE) +
    geom_point(
      data = points,
      aes(x = grid_x, y = grid_y, colour = inner_product_score_grid),
      inherit.aes = FALSE,
      size = 2.0,
      alpha = 0.95
    ) +
    scale_colour_gradient2(
      low = "#3B4CC0",
      mid = "#F7F7F7",
      high = "#B40426",
      midpoint = 0,
      limits = c(-score_limit, score_limit),
      oob = scales::squish,
      name = "Inner product score (PS)"
    ) +
    geom_segment(
      data = points,
      aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.28,
      alpha = 0.75,
      arrow = arrow(length = grid::unit(0.045, "inches"), type = "closed")
    ) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_label, y = "Dimension 2", title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.32, "cm")
    )
}

p_umap <- plot_score(cells, "umap_1", "umap_2", umap_score, "CellOracle UMAP1", "SOX4 perturbation score")
p_tsne <- plot_score(cells, "tsne_1", "tsne_2", tsne_score, "Expression t-SNE1", "SOX4 perturbation score")

save_plot <- function(plot, stem, width = 5.7, height = 4.5) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = grDevices::svg)
}
save_plot(p_umap, "figure2c_sox4_inner_product_umap")
save_plot(p_tsne, "figure2c_sox4_inner_product_tsne")
p_combined <- p_umap + p_tsne + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
save_plot(p_combined, "figure2c_sox4_inner_product_umap_tsne", width = 10.8, height = 4.5)

report <- list(
  module = "Figure 2C score-coloured companion",
  score_name = "inner_product_score_grid",
  definition = "dot(saved SOX4 perturbation flow, baseline developmental direction)",
  colour_map = "coolwarm-like blue-white-red diverging scale centered at zero",
  n_cells = nrow(cells),
  umap_grid_points = sum(umap_score$keep_score),
  tsne_grid_points = sum(tsne_score$keep_score),
  score_range_umap = range(umap_score$inner_product_score_grid[umap_score$keep_score]),
  score_range_tsne = range(tsne_score$inner_product_score_grid[tsne_score$keep_score]),
  caveat = "The project did not save a native CellOracle Gradient_calculator/Oracle_development_module grid object. These grid scores use the saved perturbation vectors and the matched R baseline direction; the t-SNE view is supplementary.",
  outputs = list(
    umap_pdf = file.path(figure_dir, "figure2c_sox4_inner_product_umap.pdf"),
    umap_png = file.path(figure_dir, "figure2c_sox4_inner_product_umap.png"),
    tsne_pdf = file.path(figure_dir, "figure2c_sox4_inner_product_tsne.pdf"),
    tsne_png = file.path(figure_dir, "figure2c_sox4_inner_product_tsne.png"),
    combined_pdf = file.path(figure_dir, "figure2c_sox4_inner_product_umap_tsne.pdf"),
    combined_png = file.path(figure_dir, "figure2c_sox4_inner_product_umap_tsne.png")
  )
)
write_json(report, file.path(data_dir, "figure2c_sox4_inner_product_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Score-coloured Figure 2C figures written to: ", figure_dir)
