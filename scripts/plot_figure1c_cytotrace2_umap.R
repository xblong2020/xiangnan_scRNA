suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1c_cytotrace2_umap.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    mode = "global",
    scores = file.path(ROOT, "metadata/figure1c/figure1c_cytotrace2_scores_by_cell.tsv.gz"),
    global_umap = file.path(ROOT, "metadata/scvi/scvi_umap.tsv.gz"),
    hep_cells = file.path(ROOT, "metadata/hepatocyte/hepatocyte_lineage_cells.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    transform = "none",
    pseudocount = 1e-4,
    max_background = 150000,
    seed = 20260708,
    width = 6.8,
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
    if (name %in% c("max_background", "seed", "width", "height", "dpi", "pseudocount")) {
      value <- as.numeric(value)
    }
    args[[name]] <- value
    i <- i + 2
  }
  args
}

read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE) && !grepl("\\.gz$", path, ignore.case = TRUE)) {
    out <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  } else {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    out <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (names(out)[1] %in% c("", "V1")) {
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

transform_scores <- function(scores, method, pseudocount) {
  raw <- as.numeric(scores)
  if (method == "none") {
    return(list(values = raw, suffix = "", label = "CytoTRACE2\nscore"))
  }
  if (method == "log10") {
    transformed <- log10(raw + pseudocount)
    rng <- range(transformed, finite = TRUE)
    if (diff(rng) <= 0) {
      scaled <- rep(0, length(transformed))
    } else {
      scaled <- (transformed - rng[1]) / diff(rng)
    }
    return(list(values = scaled, suffix = "_log10", label = "CytoTRACE2\nlog10 score"))
  }
  stop("transform must be one of: none, log10", call. = FALSE)
}

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
scores <- read_table(params$scores)

if (!"CytoTRACE2_Score" %in% names(scores)) {
  stop("Scores table must include CytoTRACE2_Score.", call. = FALSE)
}

if (params$mode == "global") {
  background <- read_table(params$global_umap)
  xy_cols <- c("UMAP_1", "UMAP_2")
  out_prefix <- "figure1C_global_cytotrace2_umap"
} else if (params$mode == "hepatocyte") {
  background <- read_table(params$hep_cells)
  xy_cols <- c("umap_hep_1", "umap_hep_2")
  out_prefix <- "figure1C_hepatocyte_lineage_cytotrace2_umap"
} else {
  stop("mode must be one of: global, hepatocyte", call. = FALSE)
}

set.seed(as.integer(params$seed))
if (is.finite(params$max_background) && nrow(background) > params$max_background) {
  background <- background[sample(seq_len(nrow(background)), size = params$max_background, replace = FALSE), , drop = FALSE]
}
background_xy <- background[, c("cell_id", xy_cols), drop = FALSE]
score_cols_keep <- setdiff(names(scores), xy_cols)
scores_clean <- scores[, score_cols_keep, drop = FALSE]

plot_df <- merge(background_xy, scores_clean, by = "cell_id", all.x = FALSE, sort = FALSE)
plot_df <- plot_df[stats::complete.cases(plot_df[, c(xy_cols, "CytoTRACE2_Score")]), , drop = FALSE]
arrows <- axis_arrows(background, xy_cols[[1]], xy_cols[[2]])
score_mapping <- transform_scores(plot_df$CytoTRACE2_Score, params$transform, params$pseudocount)
plot_df$plot_score <- score_mapping$values

palette <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")

p <- ggplot() +
  geom_point(
    data = background,
    aes(x = .data[[xy_cols[[1]]]], y = .data[[xy_cols[[2]]]]),
    color = "#D6D6D6",
    size = 0.08,
    alpha = 0.22,
    stroke = 0,
    shape = 16
  ) +
  geom_point(
    data = plot_df,
    aes(x = .data[[xy_cols[[1]]]], y = .data[[xy_cols[[2]]]], color = .data[["plot_score"]]),
    size = if (params$mode == "global") 0.18 else 0.22,
    alpha = 0.9,
    stroke = 0,
    shape = 16
  ) +
  scale_color_gradientn(colours = palette, limits = c(0, 1), name = score_mapping$label) +
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

png_path <- file.path(params$output_dir, paste0(out_prefix, score_mapping$suffix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(out_prefix, score_mapping$suffix, ".pdf"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
