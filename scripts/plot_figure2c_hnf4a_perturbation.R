#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))

root <- figure2_project_root()
data_dir <- figure2_get_arg("--data-dir", file.path(root, "metadata/driver/figure2c_hnf4a"))
b_dir <- figure2_get_arg("--figure2b-data-dir", file.path(root, "metadata/driver/figure2b_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2c_hnf4a"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")
k_neighbors <- as.integer(figure2_get_arg("--k-neighbors", "50"))
map_neighbors <- as.integer(figure2_get_arg("--map-neighbors", "40"))
seed <- as.integer(figure2_get_arg("--seed", "15071990"))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim(gzfile(file.path(data_dir, "figure2c_hnf4a_matched_cells.tsv.gz")),
                    stringsAsFactors = FALSE, check.names = FALSE)
umap_ref <- read.delim(gzfile(file.path(b_dir, "figure2b_hnf4a_baseline_grid_umap.tsv.gz")),
                       stringsAsFactors = FALSE)
tsne_ref <- read.delim(gzfile(file.path(b_dir, "figure2b_hnf4a_baseline_grid_tsne.tsv.gz")),
                       stringsAsFactors = FALSE)
cells$state_label <- factor(unname(figure2_state_labels[cells$celloracle_state]),
                            levels = unname(figure2_state_labels[figure2_state_order]))
palette_labelled <- setNames(unname(figure2_state_palette),
                             unname(figure2_state_labels[figure2_state_order]))

project_delta_to_tsne <- function(df, k = 40L) {
  k <- min(k, nrow(df) - 1L)
  out <- matrix(0, nrow(df), 2)
  for (i in seq_len(nrow(df))) {
    d <- sqrt((df$umap_1 - df$umap_1[i])^2 + (df$umap_2 - df$umap_2[i])^2)
    ix <- order(d)[seq_len(k + 1L)]
    ix <- ix[ix != i][seq_len(k)]
    dx <- df$umap_1[ix] - df$umap_1[i]
    dy <- df$umap_2[ix] - df$umap_2[i]
    bw <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2))
    design <- cbind(1, dx, dy)
    fit_x <- tryCatch(stats::lm.wfit(design, df$tsne_1[ix] - df$tsne_1[i], w = w),
                      error = function(e) NULL)
    fit_y <- tryCatch(stats::lm.wfit(design, df$tsne_2[ix] - df$tsne_2[i], w = w),
                      error = function(e) NULL)
    if (!is.null(fit_x) && !is.null(fit_y) &&
        all(is.finite(fit_x$coefficients[2:3])) && all(is.finite(fit_y$coefficients[2:3]))) {
      jacobian <- rbind(fit_x$coefficients[2:3], fit_y$coefficients[2:3])
      out[i, ] <- as.numeric(jacobian %*% c(df$delta_embedding_1[i], df$delta_embedding_2[i]))
    }
  }
  colnames(out) <- c("delta_tsne_1", "delta_tsne_2")
  as.data.frame(out)
}

cells <- cbind(cells, project_delta_to_tsne(cells, map_neighbors))
write.table(cells, gzfile(file.path(data_dir, "figure2c_hnf4a_cells_with_tsne_projection.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

aggregate_grid <- function(df, x_col, y_col, dx_col, dy_col, ref, k = 50L) {
  out <- ref
  out$flow_x <- out$flow_y <- out$flow_magnitude <- NA_real_
  for (i in seq_len(nrow(out))) {
    d <- sqrt((df[[x_col]] - out$grid_x[i])^2 + (df[[y_col]] - out$grid_y[i])^2)
    ix <- order(d)[seq_len(min(k, nrow(df)))]
    bw <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2)); w <- w / sum(w)
    out$flow_x[i] <- sum(df[[dx_col]][ix] * w)
    out$flow_y[i] <- sum(df[[dy_col]][ix] * w)
    out$flow_magnitude[i] <- sqrt(out$flow_x[i]^2 + out$flow_y[i]^2)
  }
  out$flow_norm_x <- out$flow_x / pmax(out$flow_magnitude, 1e-12)
  out$flow_norm_y <- out$flow_y / pmax(out$flow_magnitude, 1e-12)
  valid <- out$keep %in% TRUE
  cutoff <- quantile(out$flow_magnitude[valid & is.finite(out$flow_magnitude)], 0.70,
                     na.rm = TRUE, names = FALSE)
  out$show <- valid & is.finite(out$flow_magnitude) & out$flow_magnitude >= cutoff
  out$arrow_xend <- out$grid_x + out$flow_norm_x * out$arrow_length
  out$arrow_yend <- out$grid_y + out$flow_norm_y * out$arrow_length
  attr(out, "magnitude_cutoff") <- unname(cutoff)
  out
}

umap_grid <- aggregate_grid(cells, "umap_1", "umap_2", "delta_embedding_1",
                            "delta_embedding_2", umap_ref, k_neighbors)
tsne_grid <- aggregate_grid(cells, "tsne_1", "tsne_2", "delta_tsne_1",
                            "delta_tsne_2", tsne_ref, k_neighbors)
write.table(umap_grid, gzfile(file.path(data_dir, "figure2c_hnf4a_grid_umap.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne_grid, gzfile(file.path(data_dir, "figure2c_hnf4a_grid_tsne.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_field <- function(df, x_col, y_col, grid, x_lab, tag = NULL) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$show %in% TRUE, , drop = FALSE]
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = state_label)) +
    geom_point(size = 0.48, alpha = 0.30, stroke = 0) +
    geom_segment(data = arrows, aes(x = grid_x, y = grid_y, xend = arrow_xend, yend = arrow_yend),
      inherit.aes = FALSE, colour = "black", alpha = 0.82, linewidth = 0.32,
      arrow = arrow(length = grid::unit(0.05, "inches"), type = "closed")) +
    ggrepel::geom_text_repel(data = centres, aes(x = label_x, y = label_y, label = state),
      inherit.aes = FALSE, colour = "grey20", size = 2.35, family = "sans",
      min.segment.length = Inf, seed = seed) +
    scale_colour_manual(values = palette_labelled, drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = "HNF4A knockout perturbation field", tag = tag) +
    figure2_theme() + theme(legend.position = "right")
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "Figure 2C")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1")
figure2_save(p_umap, figure_dir, "figure2c_hnf4a_perturbation_umap", 5.7, 4.5)
figure2_save(p_tsne, figure_dir, "figure2c_hnf4a_perturbation_tsne", 5.7, 4.5)
figure2_save(p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right"),
             figure_dir, "figure2c_hnf4a_perturbation_umap_tsne", 10.8, 4.5)

report <- list(
  module = "Figure 2C", target_tf = target_tf,
  analysis = "HNF4A virtual knockout predicted perturbation field",
  n_cells = nrow(cells), n_grid_points = nrow(umap_grid),
  umap_arrows_shown = sum(umap_grid$show), tsne_arrows_shown = sum(tsne_grid$show),
  umap = list(method = "Direct aggregation of saved HNF4A CellOracle delta_embedding",
              k_neighbors = k_neighbors, top_magnitude_fraction = 0.30),
  tsne = list(method = "Local UMAP-to-t-SNE Jacobian projection; not native CellOracle t-SNE simulation",
              map_neighbors = map_neighbors, k_neighbors = k_neighbors),
  state_palette = as.list(figure2_state_palette),
  caveat = "Virtual knockout yields a computationally inferred state shift and network perturbation evidence."
)
figure2_write_json(report, file.path(data_dir, "figure2c_hnf4a_r_plot_report.json"))
