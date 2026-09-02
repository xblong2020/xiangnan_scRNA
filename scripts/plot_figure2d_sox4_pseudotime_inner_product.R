#!/usr/bin/env Rscript

## Figure 2D: relationship between pseudotime and SOX4 perturbation score.
## The score uses the same raw grid inner-product definition and coolwarm-like
## colour scale as the corrected Figure 2C score map.

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

c_data_dir <- normalizePath(get_arg("--figure2c-data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2c_sox4")), mustWork = FALSE)
out_dir <- normalizePath(get_arg("--out-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2d_sox4")), mustWork = FALSE)
figure_dir <- normalizePath(get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure2d_sox4")), mustWork = FALSE)
k_neighbors <- as.integer(get_arg("--k-neighbors", "50"))
n_bins <- as.integer(get_arg("--n-bins", "10"))
seed <- as.integer(get_arg("--seed", "15071990"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim(gzfile(file.path(c_data_dir, "figure2c_sox4_matched_cells.tsv.gz")), stringsAsFactors = FALSE)
umap_grid <- read.delim(gzfile(file.path(c_data_dir, "figure2c_sox4_inner_product_grid_umap.tsv.gz")), stringsAsFactors = FALSE)
tsne_grid <- read.delim(gzfile(file.path(c_data_dir, "figure2c_sox4_inner_product_grid_tsne.tsv.gz")), stringsAsFactors = FALSE)
required_cells <- c("cell_id", "pseudotime", "umap_1", "umap_2", "tsne_1", "tsne_2", "inner_product_raw_pseudotime")
missing_cells <- setdiff(required_cells, names(cells))
if (length(missing_cells)) stop("Missing matched-cell columns: ", paste(missing_cells, collapse = ", "))

## Transfer strict-main pseudotime to the same Figure 2C grid using Gaussian
## distance weights, mirroring CellOracle's grid-based comparison logic.
map_pseudotime_to_grid <- function(df, grid, x_col, y_col, k = 50L) {
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  pt <- as.numeric(df$pseudotime)
  k <- min(k, nrow(df))
  out <- grid
  out$pseudotime_grid <- NA_real_
  out$pseudotime_sd <- NA_real_
  for (i in seq_len(nrow(out))) {
    d <- sqrt((x - out$grid_x[i])^2 + (y - out$grid_y[i])^2)
    ix <- order(d)[seq_len(k)]
    bw <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2))
    w <- w / sum(w)
    mu <- sum(pt[ix] * w)
    out$pseudotime_grid[i] <- mu
    out$pseudotime_sd[i] <- sqrt(sum(w * (pt[ix] - mu)^2))
  }
  out$valid <- (out$keep_score %in% TRUE) & is.finite(out$inner_product_score_grid) & is.finite(out$pseudotime_grid)
  out
}

umap <- map_pseudotime_to_grid(cells, umap_grid, "umap_1", "umap_2", k_neighbors)
tsne <- map_pseudotime_to_grid(cells, tsne_grid, "tsne_1", "tsne_2", k_neighbors)
umap$space <- "CellOracle UMAP grid"
tsne$space <- "Expression t-SNE grid"
write.table(umap, gzfile(file.path(out_dir, "figure2d_sox4_pseudotime_inner_product_umap.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne, gzfile(file.path(out_dir, "figure2d_sox4_pseudotime_inner_product_tsne.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

summarize_bins <- function(df, space, n_bins = 10L) {
  sub <- df[df$valid, , drop = FALSE]
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  sub$pseudotime_bin <- cut(sub$pseudotime_grid, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  rows <- lapply(split(sub, sub$pseudotime_bin), function(x) {
    data.frame(
      space = space,
      pseudotime_bin = unique(x$pseudotime_bin),
      pseudotime_center = mean(x$pseudotime_grid),
      n_grid_points = nrow(x),
      score_mean = mean(x$inner_product_score_grid),
      score_median = median(x$inner_product_score_grid),
      score_q25 = unname(quantile(x$inner_product_score_grid, 0.25)),
      score_q75 = unname(quantile(x$inner_product_score_grid, 0.75)),
      positive_fraction = mean(x$inner_product_score_grid > 0)
    )
  })
  do.call(rbind, rows)
}
bin_summary <- rbind(
  summarize_bins(umap, "CellOracle UMAP grid", n_bins),
  summarize_bins(tsne, "Expression t-SNE grid", n_bins)
)
write.table(bin_summary, file.path(out_dir, "figure2d_sox4_pseudotime_bin_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Figure 2C-compatible coolwarm colours.
coolwarm_low <- "#3B4CC0"
coolwarm_mid <- "#F7F7F7"
coolwarm_high <- "#B40426"

plot_relation <- function(df, bins, title) {
  sub <- df[df$valid, , drop = FALSE]
  vm <- unname(quantile(abs(sub$inner_product_score_grid), 0.99, na.rm = TRUE))
  if (!is.finite(vm) || vm == 0) vm <- max(abs(sub$inner_product_score_grid), na.rm = TRUE)
  ggplot(sub, aes(x = pseudotime_grid, y = inner_product_score_grid)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    geom_point(aes(colour = inner_product_score_grid), size = 1.45, alpha = 0.78) +
    geom_line(
      data = bins[order(bins$pseudotime_center), ],
      aes(x = pseudotime_center, y = score_median),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.55
    ) +
    geom_point(
      data = bins,
      aes(x = pseudotime_center, y = score_median),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      stroke = 0.45,
      fill = "white",
      colour = "black"
    ) +
    scale_colour_gradient2(
      low = coolwarm_low,
      mid = coolwarm_mid,
      high = coolwarm_high,
      midpoint = 0,
      limits = c(-vm, vm),
      oob = scales::squish,
      name = "Inner product score (PS)"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Pseudotime", y = "Inner product score", title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.35, "cm")
    )
}

umap_bins <- bin_summary[bin_summary$space == "CellOracle UMAP grid", , drop = FALSE]
tsne_bins <- bin_summary[bin_summary$space == "Expression t-SNE grid", , drop = FALSE]
p_umap <- plot_relation(umap, umap_bins, "SOX4 knockout")
p_tsne <- plot_relation(tsne, tsne_bins, "SOX4 knockout (t-SNE sensitivity)")

save_plot <- function(plot, stem, width = 5.7, height = 4.3) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = grDevices::svg)
}
save_plot(p_umap, "figure2d_sox4_pseudotime_inner_product_umap")
save_plot(p_tsne, "figure2d_sox4_pseudotime_inner_product_tsne")
p_combined <- p_umap + p_tsne + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
save_plot(p_combined, "figure2d_sox4_pseudotime_inner_product_umap_tsne", width = 10.8, height = 4.3)

umap_valid <- umap[umap$valid, , drop = FALSE]
tsne_valid <- tsne[tsne$valid, , drop = FALSE]
report <- list(
  module = "Figure 2D",
  target = "SOX4 perturbation score along pseudotime",
  plotting_language = "R",
  r_version = R.version.string,
  score_definition = "raw grid inner product from corrected Figure 2C; positive means aligned with development and negative means opposing development",
  pseudotime = "driver_main_strict__pseudotime_rank mapped to grid by 50-neighbour Gaussian weighting",
  n_cells = nrow(cells),
  n_bins = n_bins,
  umap = list(
    n_grid_points = nrow(umap_valid),
    score_range = range(umap_valid$inner_product_score_grid),
    spearman_rho = unname(cor(umap_valid$pseudotime_grid, umap_valid$inner_product_score_grid, method = "spearman"))
  ),
  tsne = list(
    n_grid_points = nrow(tsne_valid),
    score_range = range(tsne_valid$inner_product_score_grid),
    spearman_rho = unname(cor(tsne_valid$pseudotime_grid, tsne_valid$inner_product_score_grid, method = "spearman"))
  ),
  colour_scale = list(low = coolwarm_low, midpoint = coolwarm_mid, high = coolwarm_high, center = 0),
  seed = seed,
  caveat = "This project-level reproduction uses the saved SOX4 perturbation vectors and matched R baseline grid. The native CellOracle Oracle_development_module object was not saved; t-SNE is a sensitivity projection.",
  outputs = list(
    main_umap_pdf = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_umap.pdf"),
    main_umap_png = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_umap.png"),
    tsne_pdf = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_tsne.pdf"),
    tsne_png = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_tsne.png"),
    combined_pdf = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_umap_tsne.pdf"),
    combined_png = file.path(figure_dir, "figure2d_sox4_pseudotime_inner_product_umap_tsne.png")
  )
)
write_json(report, file.path(out_dir, "figure2d_sox4_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Figure 2D outputs written to: ", figure_dir)
