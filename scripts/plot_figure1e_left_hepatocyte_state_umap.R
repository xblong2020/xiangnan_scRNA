suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1e_left_hepatocyte_state_umap.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    input = file.path(ROOT, "metadata/hepatocyte/hepatocyte_lineage_cells.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1E_left_hepatocyte_state_umap",
    width = 7.0,
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

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
df <- read_table(params$input)

state_order <- c(
  "normal_hepatocyte_like",
  "stressed_injured_hepatocyte",
  "regenerative_progenitor_like_hepatocyte",
  "proliferating_hepatocyte_candidate",
  "malignant_hepatocyte_candidate_needs_cnv",
  "ambiguous_epithelial_or_mixed"
)
state_labels <- c(
  normal_hepatocyte_like = "Normal hepatocyte-like",
  stressed_injured_hepatocyte = "Stressed / injured",
  regenerative_progenitor_like_hepatocyte = "Regenerative / progenitor-like",
  proliferating_hepatocyte_candidate = "Proliferating candidate",
  malignant_hepatocyte_candidate_needs_cnv = "Malignant candidate\n(needs CNV)",
  ambiguous_epithelial_or_mixed = "Ambiguous epithelial / mixed"
)
state_palette <- c(
  normal_hepatocyte_like = "#3B82A0",
  stressed_injured_hepatocyte = "#E68653",
  regenerative_progenitor_like_hepatocyte = "#49B26B",
  proliferating_hepatocyte_candidate = "#D65B9E",
  malignant_hepatocyte_candidate_needs_cnv = "#6C5CE7",
  ambiguous_epithelial_or_mixed = "#7F8C8D"
)

df$hepatocyte_state_label <- factor(df$hepatocyte_state_label, levels = state_order)
arrows <- axis_arrows(df, "umap_hep_1", "umap_hep_2")

p <- ggplot(df, aes(x = umap_hep_1, y = umap_hep_2, color = hepatocyte_state_label)) +
  geom_point(size = 0.24, alpha = 0.90, stroke = 0, shape = 16) +
  scale_color_manual(values = state_palette, labels = state_labels, name = "Hepatocyte state", drop = FALSE) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 3.6, alpha = 1, stroke = 0),
      keyheight = grid::unit(0.42, "cm"),
      keywidth = grid::unit(0.42, "cm")
    )
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
    legend.text = element_text(size = 7),
    legend.key = element_rect(fill = "white", color = NA)
  )

png_path <- file.path(params$output_dir, paste0(params$output_prefix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(params$output_prefix, ".pdf"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
