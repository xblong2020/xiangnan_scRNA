#!/usr/bin/env Rscript

## Figure 2B: SOX4 project-specific baseline developmental vector field.
## The source h5ad is prepared by prepare_figure2b_sox4_data.py; all plotting,
## t-SNE, vector-field estimation and export are performed in R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(ggrepel)
  library(Matrix)
  library(Rtsne)
  library(irlba)
  library(patchwork)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

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

data_dir <- normalizePath(
  get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2b_sox4")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure2b_sox4")),
  mustWork = FALSE
)
n_grid <- as.integer(get_arg("--n-grid", "20"))
k_neighbors <- as.integer(get_arg("--k-neighbors", "50"))
seed <- as.integer(get_arg("--seed", "15071990"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cells_path <- file.path(data_dir, "figure2b_sox4_cells.tsv.gz")
counts_path <- file.path(data_dir, "figure2b_sox4_counts.mtx.gz")
if (!file.exists(cells_path) || !file.exists(counts_path)) {
  stop("Prepared Figure 2B data are missing. Run prepare_figure2b_sox4_data.py first.")
}

cells <- read.delim(gzfile(cells_path), stringsAsFactors = FALSE, check.names = FALSE)
counts <- readMM(gzfile(counts_path))
counts <- as(counts, "dgCMatrix")
if (nrow(counts) != nrow(cells)) stop("Counts and cell metadata have different row counts")

## Calculate t-SNE from log-normalized counts, not from UMAP coordinates.
set.seed(seed)
libsize <- Matrix::rowSums(counts)
libsize[libsize <= 0] <- 1
log_counts <- Diagonal(x = 10000 / libsize) %*% counts
log_counts@x <- log1p(log_counts@x)
n_pcs <- min(50L, ncol(log_counts) - 1L, nrow(log_counts) - 1L)
pca <- irlba::prcomp_irlba(log_counts, n = n_pcs, center = TRUE, scale. = FALSE)
tsne <- Rtsne::Rtsne(
  pca$x,
  dims = 2,
  perplexity = min(30, floor((nrow(cells) - 1) / 3)),
  pca = FALSE,
  check_duplicates = FALSE,
  theta = 0.5,
  max_iter = 1000,
  verbose = FALSE,
  num_threads = 4
)$Y
cells$tsne_1 <- tsne[, 1]
cells$tsne_2 <- tsne[, 2]

write.table(
  cells,
  gzfile(file.path(data_dir, "figure2b_sox4_plot_cells.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## Estimate a local pseudotime gradient on a regular grid. The gradient is
## calculated by local weighted linear regression and normalized for plotting.
compute_grid <- function(df, x_col, y_col, value_col, n_grid = 20L, k = 50L) {
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  z <- as.numeric(df[[value_col]])
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  x <- x[ok]; y <- y[ok]; z <- z[ok]
  gx <- seq(min(x), max(x), length.out = n_grid)
  gy <- seq(min(y), max(y), length.out = n_grid)
  grid <- expand.grid(grid_x = gx, grid_y = gy)
  k <- min(k, length(x))
  kth_dist <- numeric(nrow(grid))
  beta_x <- beta_y <- rep(NA_real_, nrow(grid))
  for (i in seq_len(nrow(grid))) {
    d <- sqrt((x - grid$grid_x[i])^2 + (y - grid$grid_y[i])^2)
    ix <- order(d)[seq_len(k)]
    kth_dist[i] <- d[ix[k]]
    dx <- x[ix] - grid$grid_x[i]
    dy <- y[ix] - grid$grid_y[i]
    bandwidth <- max(stats::median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bandwidth^2))
    design <- cbind(1, dx, dy)
    fit <- tryCatch(stats::lm.wfit(design, z[ix], w = w), error = function(e) NULL)
    if (!is.null(fit) && all(is.finite(fit$coefficients[2:3]))) {
      beta_x[i] <- fit$coefficients[2]
      beta_y[i] <- fit$coefficients[3]
    }
  }
  ## Exclude low-density corners using the 70th percentile of kth-neighbor
  ## distances; this keeps arrows inside the observed cell manifold.
  cutoff <- stats::quantile(kth_dist[is.finite(kth_dist)], 0.70, names = FALSE)
  keep <- is.finite(beta_x) & is.finite(beta_y) & kth_dist <= cutoff
  norm <- sqrt(beta_x^2 + beta_y^2)
  grid$unit_x <- ifelse(keep & norm > 0, beta_x / norm, NA_real_)
  grid$unit_y <- ifelse(keep & norm > 0, beta_y / norm, NA_real_)
  grid$neighbor_distance <- kth_dist
  grid$keep <- keep
  grid$arrow_length <- max(diff(gx), diff(gy)) * 0.38
  grid$arrow_xend <- grid$grid_x + grid$unit_x * grid$arrow_length
  grid$arrow_yend <- grid$grid_y + grid$unit_y * grid$arrow_length
  attr(grid, "density_cutoff") <- unname(cutoff)
  grid
}

umap_grid <- compute_grid(cells, "umap_1", "umap_2", "pseudotime", n_grid, k_neighbors)
tsne_grid <- compute_grid(cells, "tsne_1", "tsne_2", "pseudotime", n_grid, k_neighbors)
umap_grid$space <- "CellOracle UMAP"
tsne_grid$space <- "Expression t-SNE"
write.table(umap_grid, gzfile(file.path(data_dir, "figure2b_sox4_baseline_grid_umap.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne_grid, gzfile(file.path(data_dir, "figure2b_sox4_baseline_grid_tsne.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

state_order <- c(
  "normal_reference",
  "stressed_injured",
  "regenerative_progenitor",
  "proliferating_candidate",
  "malignant_or_malignant_like"
)
state_labels <- c(
  normal_reference = "Normal/reference",
  stressed_injured = "Stressed/injured",
  regenerative_progenitor = "Regenerative/progenitor",
  proliferating_candidate = "Proliferating candidate",
  malignant_or_malignant_like = "Malignant/malignant-like"
)
cells$state_label <- factor(
  unname(state_labels[cells$celloracle_state]),
  levels = unname(state_labels[state_order])
)
pal <- ggsci::pal_lancet("lanonc")(length(state_order))
names(pal) <- unname(state_labels[state_order])

plot_field <- function(df, x_col, y_col, grid, x_lab, title) {
  centres <- aggregate(df[, c(x_col, y_col)], list(state = df$state_label), median)
  names(centres)[2:3] <- c("label_x", "label_y")
  arrows <- grid[grid$keep, , drop = FALSE]
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
    scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
    coord_equal(expand = TRUE, clip = "off") +
    labs(x = x_lab, y = "Dimension 2", title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.32, "cm"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

p_umap <- plot_field(cells, "umap_1", "umap_2", umap_grid, "CellOracle UMAP1", "Baseline developmental field")
p_tsne <- plot_field(cells, "tsne_1", "tsne_2", tsne_grid, "Expression t-SNE1", "Baseline developmental field")

save_plot <- function(plot, stem, width = 5.7, height = 4.5) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600)
  ## Use the base R SVG device so the export does not depend on svglite.
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = grDevices::svg)
}
save_plot(p_umap, "figure2b_sox4_baseline_umap")
save_plot(p_tsne, "figure2b_sox4_baseline_tsne")
p_combined <- p_umap + p_tsne + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")
save_plot(p_combined, "figure2b_sox4_baseline_umap_tsne", width = 10.8, height = 4.5)

report <- list(
  module = "Figure 2B",
  target = "SOX4 baseline developmental vector field",
  plotting_language = "R",
  r_version = R.version.string,
  palette = "ggsci::pal_lancet('lanonc')",
  n_cells = nrow(cells),
  n_genes = ncol(counts),
  pseudotime = "driver_main_strict__pseudotime_rank",
  embedding_umap = "X_celloracle_umap",
  tsne = list(method = "Rtsne on log1p library-size-normalized counts after 50-component irlba PCA", seed = seed, perplexity = min(30, floor((nrow(cells) - 1) / 3))),
  vector_field = list(method = "local weighted linear pseudotime gradient", n_grid = n_grid, k_neighbors = k_neighbors, density_quantile = 0.70, umap_density_cutoff = attr(umap_grid, "density_cutoff"), tsne_density_cutoff = attr(tsne_grid, "density_cutoff")),
  state_counts = as.list(table(cells$state_label)),
  outputs = list(
    umap_pdf = file.path(figure_dir, "figure2b_sox4_baseline_umap.pdf"),
    umap_png = file.path(figure_dir, "figure2b_sox4_baseline_umap.png"),
    umap_svg = file.path(figure_dir, "figure2b_sox4_baseline_umap.svg"),
    tsne_pdf = file.path(figure_dir, "figure2b_sox4_baseline_tsne.pdf"),
    tsne_png = file.path(figure_dir, "figure2b_sox4_baseline_tsne.png"),
    tsne_svg = file.path(figure_dir, "figure2b_sox4_baseline_tsne.svg"),
    combined_pdf = file.path(figure_dir, "figure2b_sox4_baseline_umap_tsne.pdf"),
    combined_png = file.path(figure_dir, "figure2b_sox4_baseline_umap_tsne.png"),
    combined_svg = file.path(figure_dir, "figure2b_sox4_baseline_umap_tsne.svg")
  ),
  caveat = "This is a baseline pseudotime-gradient field. SOX4 perturbation arrows belong to Figure 2C and must use the same cell subset and grid for direct comparison."
)
write_json(report, file.path(data_dir, "figure2b_sox4_r_plot_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Figure 2B UMAP and t-SNE figures written to: ", figure_dir)
