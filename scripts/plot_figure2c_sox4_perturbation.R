#!/usr/bin/env Rscript

## Figure 2C: SOX4 knockout perturbation vector field.
## UMAP uses saved CellOracle delta_embedding vectors.  The t-SNE panel is a
## supplementary local projection of the same UMAP perturbation vectors.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
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
k_neighbors <- as.integer(get_arg("--k-neighbors", "50"))
map_neighbors <- as.integer(get_arg("--map-neighbors", "40"))
seed <- as.integer(get_arg("--seed", "15071990"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cells_path <- file.path(data_dir, "figure2c_sox4_matched_cells.tsv.gz")
if (!file.exists(cells_path)) stop("Matched SOX4 data are missing. Run prepare_figure2c_sox4_data.py first.")
cells <- read.delim(gzfile(cells_path), stringsAsFactors = FALSE, check.names = FALSE)
umap_grid_path <- file.path(b_data_dir, "figure2b_sox4_baseline_grid_umap.tsv.gz")
tsne_grid_path <- file.path(b_data_dir, "figure2b_sox4_baseline_grid_tsne.tsv.gz")
if (!file.exists(umap_grid_path) || !file.exists(tsne_grid_path)) stop("Figure 2B grid data are required for matched Figure 2C plotting.")
umap_grid_ref <- read.delim(gzfile(umap_grid_path), stringsAsFactors = FALSE)
tsne_grid_ref <- read.delim(gzfile(tsne_grid_path), stringsAsFactors = FALSE)

## Use the same project palette as the existing CellOracle vector-field figures.
state_order <- c("normal_reference", "stressed_injured", "regenerative_progenitor", "proliferating_candidate", "malignant_or_malignant_like")
state_labels <- c(
  normal_reference = "Normal/reference",
  stressed_injured = "Stressed/injured",
  regenerative_progenitor = "Regenerative/progenitor",
  proliferating_candidate = "Proliferating candidate",
  malignant_or_malignant_like = "Malignant/malignant-like"
)
state_palette <- c(
  normal_reference = "#B8B8B8",
  stressed_injured = "#56B4E9",
  regenerative_progenitor = "#009E73",
  proliferating_candidate = "#E69F00",
  malignant_or_malignant_like = "#D55E00"
)
cells$state_label <- factor(unname(state_labels[cells$celloracle_state]), levels = unname(state_labels[state_order]))

## Local projection from UMAP perturbation vectors into t-SNE coordinates.
## This is required because CellOracle's saved delta_embedding is defined in
## X_celloracle_umap, while t-SNE has no native CellOracle perturbation space.
project_delta_to_tsne <- function(df, k = 40L) {
  n <- nrow(df)
  k <- min(k, n - 1L)
  out <- matrix(0, nrow = n, ncol = 2)
  for (i in seq_len(n)) {
    d <- sqrt((df$umap_1 - df$umap_1[i])^2 + (df$umap_2 - df$umap_2[i])^2)
    ix <- order(d)[seq_len(k + 1L)]
    ix <- ix[ix != i][seq_len(k)]
    dx <- df$umap_1[ix] - df$umap_1[i]
    dy <- df$umap_2[ix] - df$umap_2[i]
    bw <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2))
    design <- cbind(1, dx, dy)
    fit_x <- tryCatch(stats::lm.wfit(design, df$tsne_1[ix] - df$tsne_1[i], w = w), error = function(e) NULL)
    fit_y <- tryCatch(stats::lm.wfit(design, df$tsne_2[ix] - df$tsne_2[i], w = w), error = function(e) NULL)
    if (!is.null(fit_x) && !is.null(fit_y) && all(is.finite(fit_x$coefficients[2:3])) && all(is.finite(fit_y$coefficients[2:3]))) {
      jacobian <- rbind(fit_x$coefficients[2:3], fit_y$coefficients[2:3])
      out[i, ] <- as.numeric(jacobian %*% c(df$delta_embedding_1[i], df$delta_embedding_2[i]))
    }
  }
  colnames(out) <- c("delta_tsne_1", "delta_tsne_2")
  as.data.frame(out)
}
tsne_delta <- project_delta_to_tsne(cells, map_neighbors)
cells <- cbind(cells, tsne_delta)

aggregate_grid <- function(df, x_col, y_col, dx_col, dy_col, ref_grid, k = 50L) {
  x <- df[[x_col]]; y <- df[[y_col]]
  dx_all <- df[[dx_col]]; dy_all <- df[[dy_col]]
  k <- min(k, nrow(df))
  out <- ref_grid
  out$flow_x <- out$flow_y <- out$flow_magnitude <- NA_real_
  for (i in seq_len(nrow(out))) {
    d <- sqrt((x - out$grid_x[i])^2 + (y - out$grid_y[i])^2)
    ix <- order(d)[seq_len(k)]
    bw <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2))
    w <- w / sum(w)
    out$flow_x[i] <- sum(dx_all[ix] * w)
    out$flow_y[i] <- sum(dy_all[ix] * w)
    out$flow_magnitude[i] <- sqrt(out$flow_x[i]^2 + out$flow_y[i]^2)
  }
  out$flow_norm_x <- out$flow_x / pmax(out$flow_magnitude, 1e-12)
  out$flow_norm_y <- out$flow_y / pmax(out$flow_magnitude, 1e-12)
  out
}

umap_grid <- aggregate_grid(cells, "umap_1", "umap_2", "delta_embedding_1", "delta_embedding_2", umap_grid_ref, k_neighbors)
tsne_grid <- aggregate_grid(cells, "tsne_1", "tsne_2", "delta_tsne_1", "delta_tsne_2", tsne_grid_ref, k_neighbors)

## Match the existing CellOracle figure contract: show the strongest 30% of
## valid grid vectors, while retaining the same density mask as Figure 2B.
filter_arrows <- function(grid) {
  valid <- isTRUE(grid$keep) | grid$keep %in% TRUE
  mag <- grid$flow_magnitude
  cutoff <- stats::quantile(mag[valid & is.finite(mag)], 0.70, names = FALSE, na.rm = TRUE)
  grid$show <- valid & is.finite(mag) & mag >= cutoff
  grid$arrow_xend <- grid$grid_x + grid$flow_norm_x * grid$arrow_length
  grid$arrow_yend <- grid$grid_y + grid$flow_norm_y * grid$arrow_length
  attr(grid, "magnitude_cutoff") <- unname(cutoff)
  grid
}
umap_grid <- filter_arrows(umap_grid)
tsne_grid <- filter_arrows(tsne_grid)
write.table(umap_grid, gzfile(file.path(data_dir, "figure2c_sox4_grid_umap.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne_grid, gzfile(file.path(data_dir, "figure2c_sox4_grid_tsne.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

plot_field <- function(df, x_col, y_col, grid, x_lab, title) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$show, , drop = FALSE]
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
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = Inf,
      seed = seed
    ) +
    scale_colour_manual(values = setNames(unname(state_palette), unname(state_labels[state_order])), drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.32, "cm"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "SOX4 knockout perturbation field")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1", "SOX4 knockout perturbation field")

save_plot <- function(plot, stem, width = 5.7, height = 4.5) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = grDevices::svg)
}
save_plot(p_umap, "figure2c_sox4_perturbation_umap")
save_plot(p_tsne, "figure2c_sox4_perturbation_tsne")
p_combined <- p_umap + p_tsne + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
save_plot(p_combined, "figure2c_sox4_perturbation_umap_tsne", width = 10.8, height = 4.5)

report <- list(
  module = "Figure 2C",
  target = "SOX4 knockout perturbation vector field",
  plotting_language = "R",
  r_version = R.version.string,
  perturbation = "saved CellOracle delta_embedding after SOX4 = 0",
  n_cells = nrow(cells),
  umap_arrows_shown = sum(umap_grid$show),
  tsne_arrows_shown = sum(tsne_grid$show),
  state_palette = as.list(state_palette),
  umap = list(method = "direct aggregation of CellOracle delta_embedding", k_neighbors = k_neighbors, top_magnitude_fraction = 0.30),
  tsne = list(method = "local UMAP-to-tSNE Jacobian projection of delta_embedding", map_neighbors = map_neighbors, k_neighbors = k_neighbors, top_magnitude_fraction = 0.30),
  outputs = list(
    umap_pdf = file.path(figure_dir, "figure2c_sox4_perturbation_umap.pdf"),
    umap_png = file.path(figure_dir, "figure2c_sox4_perturbation_umap.png"),
    umap_svg = file.path(figure_dir, "figure2c_sox4_perturbation_umap.svg"),
    tsne_pdf = file.path(figure_dir, "figure2c_sox4_perturbation_tsne.pdf"),
    tsne_png = file.path(figure_dir, "figure2c_sox4_perturbation_tsne.png"),
    tsne_svg = file.path(figure_dir, "figure2c_sox4_perturbation_tsne.svg"),
    combined_pdf = file.path(figure_dir, "figure2c_sox4_perturbation_umap_tsne.pdf"),
    combined_png = file.path(figure_dir, "figure2c_sox4_perturbation_umap_tsne.png"),
    combined_svg = file.path(figure_dir, "figure2c_sox4_perturbation_umap_tsne.svg")
  ),
  caveat = "The UMAP panel is a matched 5,000-cell projection of saved CellOracle shifts. The t-SNE panel is supplementary because CellOracle perturbation vectors are defined in the CellOracle UMAP space."
)
write_json(report, file.path(data_dir, "figure2c_sox4_r_plot_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Figure 2C UMAP and t-SNE figures written to: ", figure_dir)
