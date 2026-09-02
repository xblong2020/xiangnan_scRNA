#!/usr/bin/env Rscript

## Figure 3C companion: EGR1 perturbation score on matched UMAP/t-SNE grids.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(jsonlite)
  library(scales)
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
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3c_egr1")),
  mustWork = FALSE
)
b_data_dir <- normalizePath(
  figure3_get_arg("--figure3b-data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3b_egr1")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3c_egr1")),
  mustWork = FALSE
)
k_neighbors <- as.integer(figure3_get_arg("--k-neighbors", "50"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim(
  gzfile(file.path(data_dir, "figure3c_egr1_cells_with_tsne_projection.tsv.gz")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(cells$tf == "EGR1")) stop("Figure 3C score input contains non-EGR1 rows")
umap_grid <- read.delim(gzfile(file.path(data_dir, "figure3c_egr1_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tsne_grid <- read.delim(gzfile(file.path(data_dir, "figure3c_egr1_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)
umap_ref <- read.delim(
  gzfile(file.path(b_data_dir, "figure3b_egr1_baseline_grid_umap.tsv.gz")),
  stringsAsFactors = FALSE
)
tsne_ref <- read.delim(
  gzfile(file.path(b_data_dir, "figure3b_egr1_baseline_grid_tsne.tsv.gz")),
  stringsAsFactors = FALSE
)

score_grid <- function(perturbation, reference) {
  if (nrow(perturbation) != nrow(reference) ||
      !isTRUE(all.equal(perturbation$grid_x, reference$grid_x, tolerance = 0)) ||
      !isTRUE(all.equal(perturbation$grid_y, reference$grid_y, tolerance = 0))) {
    stop("Perturbation and baseline grid coordinates are not exactly aligned")
  }
  out <- perturbation
  out$inner_product_score_grid <- (
    perturbation$flow_x * reference$unit_x +
      perturbation$flow_y * reference$unit_y
  )
  out$keep_score <- (perturbation$keep %in% TRUE) & is.finite(out$inner_product_score_grid)
  out$show_score <- out$keep_score & (perturbation$show %in% TRUE)
  out
}

compute_cell_gradient <- function(df, x_col, y_col, value_col, k = 50L) {
  coordinates <- as.matrix(df[, c(x_col, y_col)])
  values <- as.numeric(df[[value_col]])
  k <- min(k, nrow(df) - 1L)
  nn <- RANN::nn2(coordinates, query = coordinates, k = k + 1L)
  gradients <- matrix(NA_real_, nrow = nrow(df), ncol = 2)
  for (i in seq_len(nrow(df))) {
    ix_all <- nn$nn.idx[i, ]
    d_all <- nn$nn.dists[i, ]
    keep <- ix_all != i
    ix <- ix_all[keep][seq_len(k)]
    d <- d_all[keep][seq_len(k)]
    dx <- coordinates[ix, 1] - coordinates[i, 1]
    dy <- coordinates[ix, 2] - coordinates[i, 2]
    bw <- max(stats::median(d), 1e-8)
    w <- exp(-(d^2) / (2 * bw^2))
    fit <- tryCatch(stats::lm.wfit(cbind(1, dx, dy), values[ix], w = w), error = function(e) NULL)
    if (!is.null(fit) && all(is.finite(fit$coefficients[2:3]))) {
      gradients[i, ] <- fit$coefficients[2:3]
    }
  }
  norm <- sqrt(rowSums(gradients^2))
  unit <- gradients / pmax(norm, 1e-12)
  unit[!is.finite(unit)] <- NA_real_
  colnames(unit) <- c("unit_x", "unit_y")
  unit
}

umap_score <- score_grid(umap_grid, umap_ref)
tsne_score <- score_grid(tsne_grid, tsne_ref)
write.table(
  umap_score,
  gzfile(file.path(data_dir, "figure3c_egr1_inner_product_grid_umap.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  tsne_score,
  gzfile(file.path(data_dir, "figure3c_egr1_inner_product_grid_tsne.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

umap_cell_gradient <- compute_cell_gradient(cells, "umap_1", "umap_2", "pseudotime", k_neighbors)
tsne_cell_gradient <- compute_cell_gradient(cells, "tsne_1", "tsne_2", "pseudotime", k_neighbors)
cells$baseline_umap_unit_x <- umap_cell_gradient[, 1]
cells$baseline_umap_unit_y <- umap_cell_gradient[, 2]
cells$baseline_tsne_unit_x <- tsne_cell_gradient[, 1]
cells$baseline_tsne_unit_y <- tsne_cell_gradient[, 2]
cells$inner_product_score_cell_umap <- (
  cells$delta_embedding_1 * cells$baseline_umap_unit_x +
    cells$delta_embedding_2 * cells$baseline_umap_unit_y
)
cells$inner_product_score_cell_tsne <- (
  cells$delta_tsne_1 * cells$baseline_tsne_unit_x +
    cells$delta_tsne_2 * cells$baseline_tsne_unit_y
)
write.table(
  cells,
  gzfile(file.path(data_dir, "figure3c_egr1_cell_level_scores.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

score_values <- function(grid) grid$inner_product_score_grid[grid$keep_score %in% TRUE]
panel_specific_limit <- max(abs(c(score_values(umap_score), score_values(tsne_score))), na.rm = TRUE)
if (!is.finite(panel_specific_limit) || panel_specific_limit == 0) panel_specific_limit <- 1

protected_grid_paths <- c(
  file.path(PROJECT_ROOT, "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_umap.tsv.gz"),
  file.path(PROJECT_ROOT, "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_tsne.tsv.gz"),
  file.path(PROJECT_ROOT, "metadata/driver/figure2c_sox4/figure2c_sox4_inner_product_grid_umap.tsv.gz"),
  file.path(PROJECT_ROOT, "metadata/driver/figure2c_sox4/figure2c_sox4_inner_product_grid_tsne.tsv.gz")
)
shared_values <- c(score_values(umap_score), score_values(tsne_score))
shared_sources <- c()
for (path in protected_grid_paths) {
  if (!file.exists(path)) next
  table <- read.delim(gzfile(path), stringsAsFactors = FALSE)
  if (!"inner_product_score_grid" %in% names(table)) next
  valid <- if ("keep_score" %in% names(table)) table$keep_score %in% TRUE else is.finite(table$inner_product_score_grid)
  shared_values <- c(shared_values, table$inner_product_score_grid[valid & is.finite(table$inner_product_score_grid)])
  shared_sources <- c(shared_sources, figure3_norm_path(path))
}
shared_limit <- max(abs(shared_values), na.rm = TRUE)
if (!is.finite(shared_limit) || shared_limit == 0) shared_limit <- panel_specific_limit

plot_score <- function(df, x_col, y_col, grid, x_label, score_limit, tag = NULL) {
  points <- grid[grid$show_score %in% TRUE, , drop = FALSE]
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
      low = figure3_diverging[["low"]],
      mid = figure3_diverging[["mid"]],
      high = figure3_diverging[["high"]],
      midpoint = 0,
      limits = c(-score_limit, score_limit),
      oob = scales::squish,
      name = "Perturbation score (PS)"
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
    labs(x = x_label, y = "Dimension 2", title = "EGR1 perturbation score", tag = tag) +
    figure3_theme() +
    theme(legend.position = "right")
}

save_score_set <- function(limit, suffix, official = FALSE) {
  p_umap <- plot_score(
    cells, "umap_1", "umap_2", umap_score, "CellOracle UMAP1", limit,
    if (official) "Figure 3C" else NULL
  )
  p_tsne <- plot_score(cells, "tsne_1", "tsne_2", tsne_score, "Expression t-SNE1", limit)
  stem_suffix <- if (nzchar(suffix)) paste0("_", suffix) else ""
  figure3_save(p_umap, figure_dir, paste0("figure3c_egr1_inner_product_umap", stem_suffix), 5.7, 4.5)
  figure3_save(p_tsne, figure_dir, paste0("figure3c_egr1_inner_product_tsne", stem_suffix), 5.7, 4.5)
  combined <- p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right")
  figure3_save(
    combined,
    figure_dir,
    paste0("figure3c_egr1_inner_product_umap_tsne", stem_suffix),
    10.8,
    4.5
  )
}

# Main requested file names use the three-axis shared symmetric limit.
save_score_set(shared_limit, "", official = TRUE)
save_score_set(panel_specific_limit, "panel_specific")

summarize_score <- function(grid) {
  values <- score_values(grid)
  list(
    n_grid_points = length(values),
    score_range = range(values),
    positive_fraction = mean(values > 0),
    negative_fraction = mean(values < 0),
    zero_fraction = mean(values == 0)
  )
}
field_report <- jsonlite::read_json(
  file.path(data_dir, "figure3c_egr1_r_plot_report.json"),
  simplifyVector = TRUE
)
report <- list(
  module = "Figure 3C perturbation score",
  target_tf = "EGR1",
  condition = list(EGR1 = 0),
  definition = "PS = EGR1 perturbation vector dot baseline developmental vector",
  interpretation = list(
    positive = "EGR1 KO predicted displacement is aligned with the natural disease/malignant continuum.",
    negative = "EGR1 KO predicted displacement opposes the natural continuum.",
    near_zero = "No stable directional alignment was detected."
  ),
  perturbation_source = field_report$perturbation_source,
  celloracle_object = field_report$celloracle_object,
  simulation_parameters = field_report$simulation_parameters,
  propagation_steps = field_report$propagation_steps,
  n_cells = nrow(cells),
  score_summary = list(umap = summarize_score(umap_score), tsne = summarize_score(tsne_score)),
  arrow_magnitude_cutoff = field_report$arrow_magnitude_cutoff,
  colour_scale = list(
    low = unname(figure3_diverging[["low"]]),
    mid = unname(figure3_diverging[["mid"]]),
    high = unname(figure3_diverging[["high"]]),
    midpoint = 0,
    panel_specific_symmetric_limit = panel_specific_limit,
    three_axis_shared_symmetric_limit = shared_limit,
    main_files_use = "three_axis_shared_symmetric_limit"
  ),
  shared_limit_sources = as.list(shared_sources),
  methods = list(
    umap = "Direct saved EGR1 CellOracle delta_embedding aggregated on the common baseline grid.",
    tsne = "Supplementary local UMAP-to-t-SNE Jacobian projection; not a native CellOracle perturbation simulation."
  ),
  outputs = list(
    cell_level_scores = figure3_norm_path(file.path(data_dir, "figure3c_egr1_cell_level_scores.tsv.gz")),
    grid_umap = figure3_norm_path(file.path(data_dir, "figure3c_egr1_inner_product_grid_umap.tsv.gz")),
    grid_tsne = figure3_norm_path(file.path(data_dir, "figure3c_egr1_inner_product_grid_tsne.tsv.gz")),
    main_umap_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap.pdf")),
    main_umap_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap.png")),
    main_umap_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap.svg")),
    main_tsne_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_tsne.pdf")),
    main_tsne_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_tsne.png")),
    main_tsne_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_tsne.svg")),
    main_combined_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap_tsne.pdf")),
    main_combined_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap_tsne.png")),
    main_combined_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_inner_product_umap_tsne.svg"))
  ),
  caveat = "The project did not save a native CellOracle Oracle_development_module grid object. Scores use saved EGR1 displacement and the exactly matched baseline direction; t-SNE remains a supplementary projection."
)
figure3_write_json(report, file.path(data_dir, "figure3c_egr1_inner_product_report.json"))
message("Figure 3C perturbation-score outputs written to: ", figure_dir)

