suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1d_marker_feature_panel.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    plot_data = file.path(ROOT, "metadata/figure1/figure1D_marker_plot_data.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1D_marker_feature_panel",
    contour_bins = 10,
    contour_color = "#8A8A8A",
    background_color = "#D8D8D8",
    low_color = "#FFD5E4",
    high_color = "#FF2D7A",
    width = 10.5,
    height = 9.4,
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
    if (name %in% c("contour_bins", "width", "height", "dpi")) {
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

scale_gene_expression <- function(df) {
  out <- df
  out$expression_scaled <- 0
  for (gene in unique(out$gene)) {
    idx <- out$gene == gene
    vals <- as.numeric(out$expression[idx])
    if (!any(vals > 0, na.rm = TRUE)) {
      next
    }
    upper <- as.numeric(stats::quantile(vals[vals > 0], probs = 0.99, na.rm = TRUE))
    capped <- pmin(vals, upper)
    if (upper <= 0 || !is.finite(upper)) {
      scaled <- rep(0, length(vals))
    } else {
      scaled <- capped / upper
    }
    out$expression_scaled[idx] <- scaled
  }
  out
}

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_df <- read_table(params$plot_data)
plot_df <- scale_gene_expression(plot_df)
plot_df$gene <- factor(plot_df$gene, levels = unique(plot_df$gene))

background_df <- unique(plot_df[, c("cell_id", "UMAP_1", "UMAP_2", "gene"), drop = FALSE])
expr_df <- plot_df[plot_df$expression_scaled > 0, , drop = FALSE]

p <- ggplot() +
  geom_point(
    data = background_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = params$background_color,
    size = 0.08,
    alpha = 0.24,
    stroke = 0,
    shape = 16
  ) +
  stat_density_2d(
    data = background_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = params$contour_color,
    linewidth = 0.20,
    bins = as.integer(params$contour_bins),
    alpha = 0.75
  ) +
  geom_point(
    data = expr_df,
    aes(x = UMAP_1, y = UMAP_2, color = expression_scaled),
    size = 0.12,
    alpha = 0.95,
    stroke = 0,
    shape = 16
  ) +
  scale_color_gradient(low = params$low_color, high = params$high_color, limits = c(0, 1), name = "Scaled\nexpression") +
  facet_wrap(~gene, ncol = 3) +
  coord_equal() +
  theme_void(base_size = 8) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.spacing = grid::unit(0.55, "lines")
  )

png_path <- file.path(params$output_dir, paste0(params$output_prefix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(params$output_prefix, ".pdf"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
