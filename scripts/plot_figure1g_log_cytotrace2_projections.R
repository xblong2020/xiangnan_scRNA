suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1g_log_cytotrace2_projections.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
source(file.path(ROOT, "scripts", "figure1g_stemness_helpers.R"))

parse_args <- function(argv) {
  defaults <- list(
    scores = file.path(ROOT, "metadata/figure1c/figure1c_cytotrace2_scores_by_cell.hepatocyte.tsv.gz"),
    tsne = file.path(ROOT, "metadata/figure1/figure1g_hepatocyte_cytotrace2_tsne.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1G_hepatocyte_lineage_cytotrace2_log10",
    pseudocount = 1e-4,
    width = 6.8,
    height = 6.2,
    dpi = 600
  )
  args <- defaults
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key, call. = FALSE)
    name <- sub("^--", "", key)
    if (!name %in% names(args)) stop("Unknown argument: --", name, call. = FALSE)
    if (i == length(argv)) stop("Missing value for --", name, call. = FALSE)
    value <- argv[[i + 1L]]
    if (name %in% c("pseudocount", "width", "height", "dpi")) value <- as.numeric(value)
    args[[name]] <- value
    i <- i + 2L
  }
  args
}

read_table <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (names(out)[[1]] %in% c("", "V1", "Unnamed: 0")) names(out)[[1]] <- "cell_id"
  out
}

axis_arrows <- function(df, x_col, y_col, x_label, y_label) {
  x_range <- range(df[[x_col]], finite = TRUE)
  y_range <- range(df[[y_col]], finite = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  x0 <- x_range[[1]] + 0.05 * x_span
  y0 <- y_range[[1]] + 0.06 * y_span
  data.frame(
    axis = c(x_label, y_label),
    x = c(x0, x0), y = c(y0, y0),
    xend = c(x0 + 0.12 * x_span, x0), yend = c(y0, y0 + 0.12 * y_span),
    label_x = c(x0 + 0.14 * x_span, x0 - 0.02 * x_span),
    label_y = c(y0 - 0.015 * y_span, y0 + 0.14 * y_span),
    stringsAsFactors = FALSE
  )
}

build_projection_plot <- function(background, scored, x_col, y_col, title, x_label, y_label) {
  arrows <- axis_arrows(background, x_col, y_col, x_label, y_label)
  palette <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  ggplot() +
    geom_point(
      data = background,
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      color = "#D6D6D6", size = 0.16, alpha = 0.40, shape = 16, stroke = 0
    ) +
    geom_density_2d(
      data = background,
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      color = "#A7A7A7", linewidth = 0.25, bins = 8, alpha = 0.65
    ) +
    geom_point(
      data = scored,
      aes(x = .data[[x_col]], y = .data[[y_col]], color = log10_cytotrace2_scaled),
      size = 0.48, alpha = 0.90, shape = 16, stroke = 0
    ) +
    scale_color_gradientn(
      colours = palette,
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("low", "mid", "high"),
      name = "Log10-normalized CytoTRACE2 score"
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = "black", linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(0.09, "in"), type = "closed")
    ) +
    geom_text(
      data = arrows,
      aes(x = label_x, y = label_y, label = axis),
      inherit.aes = FALSE, color = "black", size = 2.4
    ) +
    labs(title = title) +
    coord_equal(clip = "off") +
    theme_void(base_size = 8) +
    theme(
      plot.margin = margin(6, 8, 6, 6),
      plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 3)),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
score_source <- read_table(params$scores)
scores <- prepare_figure1g_data(score_source, unique(score_source$hepatocyte_state_label))
required_score_umap <- c("cell_id", "umap_hep_1", "umap_hep_2")
if (length(setdiff(required_score_umap, names(score_source))) > 0L) {
  stop("CytoTRACE2 score table must include cell_id, umap_hep_1, and umap_hep_2.", call. = FALSE)
}
scores <- merge(scores, score_source[, required_score_umap], by = "cell_id", all.x = TRUE, sort = FALSE)
scores$log10_cytotrace2_scaled <- log10_minmax_scores(scores$CytoTRACE2_Score, params$pseudocount)
umap <- scores[, c("cell_id", "umap_hep_1", "umap_hep_2"), drop = FALSE]
names(umap)[names(umap) == "umap_hep_1"] <- "UMAP_1"
names(umap)[names(umap) == "umap_hep_2"] <- "UMAP_2"
tsne <- read_table(params$tsne)

required_umap <- c("cell_id", "UMAP_1", "UMAP_2")
required_tsne <- c("cell_id", "TSNE_1", "TSNE_2")
if (length(setdiff(required_umap, names(umap))) > 0L) stop("UMAP table must include cell_id, UMAP_1, and UMAP_2.", call. = FALSE)
if (length(setdiff(required_tsne, names(tsne))) > 0L) stop("t-SNE table must include cell_id, TSNE_1, and TSNE_2.", call. = FALSE)

umap_scored <- merge(umap, scores[, c("cell_id", "CytoTRACE2_Score", "log10_cytotrace2_scaled")], by = "cell_id", all = FALSE, sort = FALSE)
tsne_scored <- merge(tsne, scores[, c("cell_id", "CytoTRACE2_Score", "log10_cytotrace2_scaled")], by = "cell_id", all = FALSE, sort = FALSE)
if (nrow(umap_scored) == 0L || nrow(tsne_scored) == 0L) stop("No cells overlap between trajectory coordinates and CytoTRACE2 scores.", call. = FALSE)

umap_plot <- build_projection_plot(umap, umap_scored, "UMAP_1", "UMAP_2", "Figure 1G. CytoTRACE2 stemness projection (UMAP)", "UMAP1", "UMAP2")
tsne_plot <- build_projection_plot(tsne, tsne_scored, "TSNE_1", "TSNE_2", "Figure 1G. CytoTRACE2 stemness projection (t-SNE)", "t-SNE1", "t-SNE2")

umap_png <- file.path(params$output_dir, paste0(params$output_prefix, "_umap.png"))
umap_pdf <- file.path(params$output_dir, paste0(params$output_prefix, "_umap.pdf"))
tsne_png <- file.path(params$output_dir, paste0(params$output_prefix, "_tsne.png"))
tsne_pdf <- file.path(params$output_dir, paste0(params$output_prefix, "_tsne.pdf"))
data_path <- file.path(params$output_dir, paste0(params$output_prefix, "_plot_data.tsv.gz"))
ggsave(umap_png, umap_plot, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(umap_pdf, umap_plot, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
ggsave(tsne_png, tsne_plot, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(tsne_pdf, tsne_plot, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
plot_data <- merge(umap_scored, tsne[, required_tsne], by = "cell_id", all = FALSE, sort = FALSE)
con <- gzfile(data_path, open = "wt")
on.exit(close(con), add = TRUE)
write.table(plot_data, con, sep = "\t", quote = FALSE, row.names = FALSE)
message("WROTE ", normalizePath(umap_png, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(tsne_png, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(data_path, winslash = "/", mustWork = FALSE))
