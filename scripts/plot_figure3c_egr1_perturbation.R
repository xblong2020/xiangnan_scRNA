#!/usr/bin/env Rscript

## Figure 3C: EGR1 knockout perturbation vector field.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(RANN)
  library(jsonlite)
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
map_neighbors <- as.integer(figure3_get_arg("--map-neighbors", "40"))
seed <- as.integer(figure3_get_arg("--seed", "15071990"))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cells_path <- file.path(data_dir, "figure3c_egr1_matched_cells.tsv.gz")
umap_grid_path <- file.path(b_data_dir, "figure3b_egr1_baseline_grid_umap.tsv.gz")
tsne_grid_path <- file.path(b_data_dir, "figure3b_egr1_baseline_grid_tsne.tsv.gz")
if (!all(file.exists(c(cells_path, umap_grid_path, tsne_grid_path)))) {
  stop("Figure 3B grids and matched EGR1 cell vectors are required.")
}
cells <- read.delim(gzfile(cells_path), stringsAsFactors = FALSE, check.names = FALSE)
if (!all(cells$tf == "EGR1")) stop("Figure 3C input contains non-EGR1 perturbation rows")
umap_grid_ref <- read.delim(gzfile(umap_grid_path), stringsAsFactors = FALSE)
tsne_grid_ref <- read.delim(gzfile(tsne_grid_path), stringsAsFactors = FALSE)
cells$state_label <- factor(
  unname(figure3_state_labels[cells$celloracle_state]),
  levels = unname(figure3_state_labels[figure3_state_order])
)
palette_labelled <- setNames(
  unname(figure3_state_palette),
  unname(figure3_state_labels[figure3_state_order])
)

project_delta_to_tsne <- function(df, k = 40L) {
  umap <- as.matrix(df[, c("umap_1", "umap_2")])
  tsne <- as.matrix(df[, c("tsne_1", "tsne_2")])
  k <- min(k, nrow(df) - 1L)
  nn <- RANN::nn2(umap, query = umap, k = k + 1L)
  out <- matrix(0, nrow = nrow(df), ncol = 2)
  for (i in seq_len(nrow(df))) {
    ix_all <- nn$nn.idx[i, ]
    dist_all <- nn$nn.dists[i, ]
    keep <- ix_all != i
    ix <- ix_all[keep][seq_len(k)]
    d <- dist_all[keep][seq_len(k)]
    dx <- umap[ix, 1] - umap[i, 1]
    dy <- umap[ix, 2] - umap[i, 2]
    bw <- max(stats::median(d), 1e-8)
    w <- exp(-(d^2) / (2 * bw^2))
    design <- cbind(1, dx, dy)
    fit_x <- tryCatch(stats::lm.wfit(design, tsne[ix, 1] - tsne[i, 1], w = w), error = function(e) NULL)
    fit_y <- tryCatch(stats::lm.wfit(design, tsne[ix, 2] - tsne[i, 2], w = w), error = function(e) NULL)
    if (!is.null(fit_x) && !is.null(fit_y) &&
        all(is.finite(fit_x$coefficients[2:3])) && all(is.finite(fit_y$coefficients[2:3]))) {
      jacobian <- rbind(fit_x$coefficients[2:3], fit_y$coefficients[2:3])
      out[i, ] <- as.numeric(jacobian %*% c(df$delta_embedding_1[i], df$delta_embedding_2[i]))
    }
  }
  colnames(out) <- c("delta_tsne_1", "delta_tsne_2")
  as.data.frame(out)
}

tsne_delta <- project_delta_to_tsne(cells, map_neighbors)
cells <- cbind(cells, tsne_delta)
write.table(
  cells,
  gzfile(file.path(data_dir, "figure3c_egr1_cells_with_tsne_projection.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

aggregate_grid <- function(df, x_col, y_col, dx_col, dy_col, ref_grid, k = 50L) {
  coordinates <- as.matrix(df[, c(x_col, y_col)])
  query <- as.matrix(ref_grid[, c("grid_x", "grid_y")])
  k <- min(k, nrow(df))
  nn <- RANN::nn2(coordinates, query = query, k = k)
  out <- ref_grid
  out$flow_x <- out$flow_y <- out$flow_magnitude <- NA_real_
  for (i in seq_len(nrow(out))) {
    ix <- nn$nn.idx[i, ]
    d <- nn$nn.dists[i, ]
    bw <- max(stats::median(d), 1e-8)
    w <- exp(-(d^2) / (2 * bw^2))
    w <- w / sum(w)
    out$flow_x[i] <- sum(df[[dx_col]][ix] * w)
    out$flow_y[i] <- sum(df[[dy_col]][ix] * w)
    out$flow_magnitude[i] <- sqrt(out$flow_x[i]^2 + out$flow_y[i]^2)
  }
  out$flow_norm_x <- out$flow_x / pmax(out$flow_magnitude, 1e-12)
  out$flow_norm_y <- out$flow_y / pmax(out$flow_magnitude, 1e-12)
  out
}

filter_arrows <- function(grid) {
  valid <- grid$keep %in% TRUE
  magnitude <- grid$flow_magnitude
  cutoff <- stats::quantile(magnitude[valid & is.finite(magnitude)], 0.70, names = FALSE, na.rm = TRUE)
  grid$show <- valid & is.finite(magnitude) & magnitude >= cutoff
  grid$arrow_xend <- grid$grid_x + grid$flow_norm_x * grid$arrow_length
  grid$arrow_yend <- grid$grid_y + grid$flow_norm_y * grid$arrow_length
  attr(grid, "magnitude_cutoff") <- unname(cutoff)
  grid
}

umap_grid <- filter_arrows(
  aggregate_grid(cells, "umap_1", "umap_2", "delta_embedding_1", "delta_embedding_2", umap_grid_ref, k_neighbors)
)
tsne_grid <- filter_arrows(
  aggregate_grid(cells, "tsne_1", "tsne_2", "delta_tsne_1", "delta_tsne_2", tsne_grid_ref, k_neighbors)
)
write.table(
  umap_grid,
  gzfile(file.path(data_dir, "figure3c_egr1_grid_umap.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  tsne_grid,
  gzfile(file.path(data_dir, "figure3c_egr1_grid_tsne.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_field <- function(df, x_col, y_col, grid, x_lab, tag = NULL) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$show %in% TRUE, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = state_label)) +
    geom_point(size = 0.48, alpha = 0.30, stroke = 0) +
    geom_segment(
      data = arrows,
      aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE,
      colour = "black",
      alpha = 0.82,
      linewidth = 0.32,
      arrow = arrow(length = grid::unit(0.05, "inches"), type = "closed")
    ) +
    ggrepel::geom_text_repel(
      data = centres,
      aes(x = label_x, y = label_y, label = state),
      inherit.aes = FALSE,
      colour = "grey20",
      size = 2.35,
      family = "sans",
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = Inf,
      seed = seed
    ) +
    scale_colour_manual(values = palette_labelled, drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = "EGR1 knockout perturbation field", tag = tag) +
    figure3_theme() +
    theme(legend.position = "right")
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "Figure 3C")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1")
figure3_save(p_umap, figure_dir, "figure3c_egr1_perturbation_umap", 5.7, 4.5)
figure3_save(p_tsne, figure_dir, "figure3c_egr1_perturbation_tsne", 5.7, 4.5)
p_combined <- p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right")
figure3_save(p_combined, figure_dir, "figure3c_egr1_perturbation_umap_tsne", 10.8, 4.5)

data_report_path <- file.path(data_dir, "figure3c_egr1_data_report.json")
data_report <- jsonlite::read_json(data_report_path, simplifyVector = TRUE)
report <- list(
  module = "Figure 3C perturbation field",
  target_tf = "EGR1",
  condition = list(EGR1 = 0),
  perturbation_source = data_report$source_perturbation,
  celloracle_object = data_report$source_oracle,
  simulation_parameters = data_report$simulation_parameters,
  propagation_steps = data_report$simulation_parameters$n_propagation,
  n_cells = nrow(cells),
  state_palette = as.list(figure3_state_palette),
  arrow_magnitude_cutoff = list(
    umap = attr(umap_grid, "magnitude_cutoff"),
    tsne = attr(tsne_grid, "magnitude_cutoff")
  ),
  arrows_shown = list(umap = sum(umap_grid$show), tsne = sum(tsne_grid$show)),
  umap = list(
    method = "Direct Gaussian aggregation of saved EGR1 CellOracle delta_embedding on the exact Figure 3B grid",
    k_neighbors = k_neighbors,
    top_magnitude_fraction = 0.30
  ),
  tsne = list(
    method = "Local UMAP-to-t-SNE Jacobian projection followed by Gaussian aggregation",
    map_neighbors = map_neighbors,
    k_neighbors = k_neighbors,
    top_magnitude_fraction = 0.30,
    native_celloracle_simulation = FALSE
  ),
  outputs = list(
    cell_level = figure3_norm_path(file.path(data_dir, "figure3c_egr1_cells_with_tsne_projection.tsv.gz")),
    grid_umap = figure3_norm_path(file.path(data_dir, "figure3c_egr1_grid_umap.tsv.gz")),
    grid_tsne = figure3_norm_path(file.path(data_dir, "figure3c_egr1_grid_tsne.tsv.gz")),
    umap_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap.pdf")),
    umap_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap.png")),
    umap_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap.svg")),
    tsne_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_tsne.pdf")),
    tsne_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_tsne.png")),
    tsne_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_tsne.svg")),
    combined_pdf = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap_tsne.pdf")),
    combined_png = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap_tsne.png")),
    combined_svg = figure3_norm_path(file.path(figure_dir, "figure3c_egr1_perturbation_umap_tsne.svg"))
  ),
  caveat = "t-SNE is a supplementary projection and not a native CellOracle perturbation simulation."
)
figure3_write_json(report, file.path(data_dir, "figure3c_egr1_r_plot_report.json"))
message("Figure 3C perturbation field written to: ", figure_dir)

