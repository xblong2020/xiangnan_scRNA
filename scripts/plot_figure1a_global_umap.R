suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1a_global_umap.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    umap = file.path(ROOT, "metadata/scvi/scvi_umap.tsv.gz"),
    annotations = file.path(ROOT, "metadata/celltype/celltypist_major_by_cell.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1A_global_integrated_umap",
    label_col = "major_celltype",
    target_lineage = "",
    point_size = 0.12,
    point_alpha = 0.85,
    seed = 20260601,
    max_cells = Inf,
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
    if (name %in% c("point_size", "point_alpha", "width", "height", "dpi", "max_cells", "seed")) {
      value <- as.numeric(value)
    }
    args[[name]] <- value
    i <- i + 2
  }
  args
}

read_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  if (requireNamespace("data.table", quietly = TRUE) && !grepl("\\.gz$", path, ignore.case = TRUE)) {
    out <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  } else {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    out <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  first_name <- names(out)[[1]]
  if (first_name == "" || first_name == "V1") {
    names(out)[[1]] <- "cell_id"
  }
  out
}

make_lancet_palette <- function(labels) {
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop("R package ggsci is required for the Lancet palette.", call. = FALSE)
  }
  base <- ggsci::pal_lancet("lanonc")(9)
  cols <- if (length(labels) <= length(base)) base[seq_along(labels)] else grDevices::colorRampPalette(base)(length(labels))
  stats::setNames(cols, labels)
}

convex_hull <- function(df, label_col, target_lineage) {
  if (is.na(target_lineage) || target_lineage == "" || !target_lineage %in% df[[label_col]]) {
    return(NULL)
  }
  target <- df[df[[label_col]] == target_lineage, c("UMAP_1", "UMAP_2"), drop = FALSE]
  target <- target[stats::complete.cases(target), , drop = FALSE]
  if (nrow(target) < 3) {
    return(NULL)
  }
  target[grDevices::chull(target$UMAP_1, target$UMAP_2), , drop = FALSE]
}

axis_arrows <- function(df) {
  x_rng <- range(df$UMAP_1, finite = TRUE)
  y_rng <- range(df$UMAP_2, finite = TRUE)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(as.integer(args$seed))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

umap <- read_table(args$umap)
annot <- read_table(args$annotations)

required_umap <- c("cell_id", "UMAP_1", "UMAP_2")
missing_umap <- setdiff(required_umap, names(umap))
if (length(missing_umap) > 0) {
  stop("UMAP table is missing required columns: ", paste(missing_umap, collapse = ", "), call. = FALSE)
}
required_annot <- c("cell_id", args$label_col)
missing_annot <- setdiff(required_annot, names(annot))
if (length(missing_annot) > 0) {
  stop("Annotation table is missing required columns: ", paste(missing_annot, collapse = ", "), call. = FALSE)
}

plot_df <- merge(
  umap[, required_umap, drop = FALSE],
  annot[, required_annot, drop = FALSE],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
plot_df[[args$label_col]][is.na(plot_df[[args$label_col]]) | plot_df[[args$label_col]] == ""] <- "Unknown"
plot_df <- plot_df[stats::complete.cases(plot_df[, c("UMAP_1", "UMAP_2")]), , drop = FALSE]

if (is.finite(args$max_cells) && nrow(plot_df) > args$max_cells) {
  keep <- sample(seq_len(nrow(plot_df)), size = args$max_cells, replace = FALSE)
  plot_df <- plot_df[sort(keep), , drop = FALSE]
}

counts <- sort(table(plot_df[[args$label_col]]), decreasing = TRUE)
label_levels <- names(counts)
plot_df[[args$label_col]] <- factor(plot_df[[args$label_col]], levels = label_levels)
palette <- make_lancet_palette(label_levels)
hull <- convex_hull(plot_df, args$label_col, args$target_lineage)
arrows <- axis_arrows(plot_df)

p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[args$label_col]])) +
  geom_point(size = args$point_size, alpha = args$point_alpha, stroke = 0, shape = 16) +
  scale_color_manual(values = palette, name = "Cell Type") +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), keyheight = grid::unit(0.34, "cm"))) +
  coord_equal(clip = "off") +
  theme_void(base_size = 8) +
  theme(
    plot.margin = margin(6, 6, 6, 6),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.4),
    legend.key = element_blank()
  )

if (!is.null(hull)) {
  p <- p +
    geom_polygon(
      data = hull,
      aes(x = UMAP_1, y = UMAP_2),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.35,
      linetype = "22"
    )
}

p <- p +
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
  )

png_path <- file.path(args$output_dir, paste0(args$output_prefix, ".png"))
pdf_path <- file.path(args$output_dir, paste0(args$output_prefix, ".pdf"))
ggsave(png_path, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
ggsave(pdf_path, p, width = args$width, height = args$height, bg = "white", useDingbats = FALSE)

counts_path <- file.path(args$output_dir, paste0(args$output_prefix, "_celltype_counts.tsv"))
write.table(
  data.frame(cell_type = names(counts), n_cells = as.integer(counts), row.names = NULL),
  counts_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(counts_path, winslash = "/", mustWork = FALSE))
