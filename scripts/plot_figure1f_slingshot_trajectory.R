suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1f_slingshot_trajectory.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    embedding = "umap",
    pseudotime = file.path(ROOT, "metadata/trajectory/trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz"),
    cell_metadata = file.path(ROOT, "data/processed/trajectory/module5_3/main_strict/cell_metadata.tsv.gz"),
    umap = file.path(ROOT, "data/processed/trajectory/module5_3/main_strict/embedding_umap.tsv.gz"),
    tsne = file.path(ROOT, "metadata/trajectory/figure1f_main_strict_tsne.tsv.gz"),
    cluster_summary = file.path(ROOT, "data/processed/trajectory/module5_3/main_strict/slingshot_cluster_summary.tsv"),
    output_dir = file.path(ROOT, "figures/figure1"),
    width = 7.2,
    height = 6.2,
    dpi = 600
  )
  args <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(args)) {
      stop("Unknown argument: --", name, call. = FALSE)
    }
    if (i == length(argv)) {
      stop("Missing value for --", name, call. = FALSE)
    }
    value <- argv[[i + 1]]
    if (name %in% c("width", "height", "dpi")) {
      value <- as.numeric(value)
    }
    args[[name]] <- value
    i <- i + 2
  }
  args
}

read_table <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (names(out)[1] %in% c("", "V1", "Unnamed: 0")) {
    names(out)[1] <- "cell_id"
  }
  out
}

axis_arrows <- function(df, x_col, y_col) {
  x_rng <- range(df[[x_col]], finite = TRUE)
  y_rng <- range(df[[y_col]], finite = TRUE)
  x_span <- diff(x_rng)
  y_span <- diff(y_rng)
  x0 <- x_rng[[1]] + 0.05 * x_span
  y0 <- y_rng[[1]] + 0.06 * y_span
  data.frame(
    axis = c("UMAP1", "UMAP2"),
    x = c(x0, x0),
    y = c(y0, y0),
    xend = c(x0 + 0.12 * x_span, x0),
    yend = c(y0, y0 + 0.12 * y_span),
    label_x = c(x0 + 0.14 * x_span, x0 - 0.02 * x_span),
    label_y = c(y0 - 0.015 * y_span, y0 + 0.14 * y_span),
    stringsAsFactors = FALSE
  )
}

fit_branch_curve <- function(df, x_col, y_col, pt_col) {
  df <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]) & is.finite(df[[pt_col]]), , drop = FALSE]
  df <- df[order(df[[pt_col]]), , drop = FALSE]
  df <- df[!duplicated(df[[pt_col]]), , drop = FALSE]
  if (nrow(df) < 4) {
    return(NULL)
  }
  grid <- seq(min(df[[pt_col]]), max(df[[pt_col]]), length.out = 200)
  sp_x <- stats::smooth.spline(x = df[[pt_col]], y = df[[x_col]], spar = 0.55)
  sp_y <- stats::smooth.spline(x = df[[pt_col]], y = df[[y_col]], spar = 0.55)
  out <- data.frame(
    pseudotime = grid,
    x = stats::predict(sp_x, x = grid)$y,
    y = stats::predict(sp_y, x = grid)$y
  )
  out[is.finite(out$x) & is.finite(out$y), , drop = FALSE]
}

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)

pt <- read_table(params$pseudotime)
meta <- read_table(params$cell_metadata)
cluster_summary <- read_table(params$cluster_summary)

if (params$embedding == "umap") {
  embed <- read_table(params$umap)
  x_col <- "UMAP_1"
  y_col <- "UMAP_2"
  out_prefix <- "figure1F_slingshot_trajectory_umap"
} else if (params$embedding == "tsne") {
  embed <- read_table(params$tsne)
  x_col <- "TSNE_1"
  y_col <- "TSNE_2"
  out_prefix <- "figure1F_slingshot_trajectory_tsne"
} else {
  stop("embedding must be one of: umap, tsne", call. = FALSE)
}

meta_keep <- meta[, setdiff(names(meta), c("slingshot_cluster", "trajectory_root_end_role", "cell_disease_stage", "sample_disease_stage")), drop = FALSE]
plot_df <- merge(pt, meta_keep, by = "cell_id", all.x = TRUE, sort = FALSE)
plot_df <- merge(plot_df, embed, by = "cell_id", all.x = TRUE, sort = FALSE)
pt_col <- "main_strict__slingshot_scanvi_norm"
if (!pt_col %in% names(plot_df)) {
  stop("Missing slingshot pseudotime column: ", pt_col, call. = FALSE)
}

cluster_medians <- aggregate(
  cbind(main_strict__slingshot_scanvi_norm, main_strict__slingshot_hepatocyte_pca_norm) ~ slingshot_cluster,
  data = plot_df,
  FUN = median,
  na.rm = TRUE
)
cluster_summary$cluster <- as.character(cluster_summary$cluster)
cluster_medians$slingshot_cluster <- as.character(cluster_medians$slingshot_cluster)
cluster_summary <- merge(cluster_summary, cluster_medians, by.x = "cluster", by.y = "slingshot_cluster", all.x = TRUE, sort = FALSE)

root_clusters <- cluster_summary$cluster[cluster_summary$root_fraction >= 0.5]
trunk_clusters <- cluster_summary$cluster[
  cluster_summary$root_fraction >= 0.5 |
    (
      cluster_summary$malignant_fraction < 0.1 &
        cluster_summary$progenitor_fraction < 0.1 &
        cluster_summary$main_strict__slingshot_scanvi_norm <= 0.53
    )
]
malignant_clusters <- unique(c(trunk_clusters, cluster_summary$cluster[cluster_summary$malignant_fraction >= 0.5]))
progenitor_clusters <- unique(c(trunk_clusters, cluster_summary$cluster[cluster_summary$progenitor_fraction >= 0.5]))

plot_df$slingshot_cluster <- as.character(plot_df$slingshot_cluster)
cluster_positions <- aggregate(
  cbind(x_coord = plot_df[[x_col]], y_coord = plot_df[[y_col]], pseudotime = plot_df[[pt_col]]) ~ slingshot_cluster,
  data = plot_df,
  FUN = median,
  na.rm = TRUE
)
malignant_df <- cluster_positions[cluster_positions$slingshot_cluster %in% malignant_clusters, , drop = FALSE]
progenitor_df <- cluster_positions[cluster_positions$slingshot_cluster %in% progenitor_clusters, , drop = FALSE]
malignant_curve <- fit_branch_curve(malignant_df, "x_coord", "y_coord", "pseudotime")
progenitor_curve <- fit_branch_curve(progenitor_df, "x_coord", "y_coord", "pseudotime")
if (!is.null(malignant_curve)) {
  malignant_curve$branch <- "Malignant branch"
}
if (!is.null(progenitor_curve)) {
  progenitor_curve$branch <- "Progenitor branch"
}
curve_df <- rbind(malignant_curve, progenitor_curve)

plot_df <- plot_df[is.finite(plot_df[[x_col]]) & is.finite(plot_df[[y_col]]) & is.finite(plot_df[[pt_col]]), , drop = FALSE]
arrows <- axis_arrows(plot_df, x_col, y_col)
branch_colors <- c("Malignant branch" = "#F28E2B", "Progenitor branch" = "#4E79A7")
pt_palette <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")

p <- ggplot() +
  geom_point(
    data = plot_df,
    aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[pt_col]]),
    shape = 16,
    size = 0.45,
    alpha = 0.82,
    stroke = 0
  ) +
  scale_color_gradientn(colours = pt_palette, limits = c(0, 1), name = "Pseudotime") +
  geom_path(
    data = curve_df[curve_df$branch == "Malignant branch", , drop = FALSE],
    aes(x = x, y = y, group = branch),
    linewidth = 1.15,
    inherit.aes = FALSE,
    color = "#F28E2B"
  ) +
  geom_path(
    data = curve_df[curve_df$branch == "Progenitor branch", , drop = FALSE],
    aes(x = x, y = y, group = branch),
    linewidth = 1.15,
    inherit.aes = FALSE,
    color = "#4E79A7"
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35,
    arrow = grid::arrow(length = grid::unit(0.09, "in"), type = "closed")
  ) +
  geom_text(
    data = arrows,
    aes(x = label_x, y = label_y, label = axis),
    inherit.aes = FALSE,
    color = "black",
    size = 2.4
  ) +
  coord_equal(clip = "off") +
  theme_void(base_size = 8) +
  theme(
    plot.margin = margin(6, 6, 6, 6),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

if (!is.null(malignant_curve) && nrow(malignant_curve) > 0) {
  malignant_end <- malignant_curve[nrow(malignant_curve), , drop = FALSE]
  p <- p + annotate("text", x = malignant_end$x, y = malignant_end$y, label = "Malignant", color = "#F28E2B", size = 2.8, hjust = -0.05)
}
if (!is.null(progenitor_curve) && nrow(progenitor_curve) > 0) {
  progenitor_end <- progenitor_curve[nrow(progenitor_curve), , drop = FALSE]
  p <- p + annotate("text", x = progenitor_end$x, y = progenitor_end$y, label = "Progenitor", color = "#4E79A7", size = 2.8, hjust = -0.05)
}

png_path <- file.path(params$output_dir, paste0(out_prefix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(out_prefix, ".pdf"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
