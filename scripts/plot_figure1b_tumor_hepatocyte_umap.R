suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1b_tumor_hepatocyte_umap.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    umap = file.path(ROOT, "metadata/scvi/scvi_umap.tsv.gz"),
    annotations = file.path(ROOT, "metadata/celltype/celltypist_major_by_cell.tsv.gz"),
    sample_summary = file.path(ROOT, "metadata/sample_info/sample_summary/all_datasets_sample_summary.tsv"),
    plot_data = file.path(ROOT, "metadata/figure1/figure1B_tumor_hepatocyte_plot_data.tsv"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1B_tumor_hepatocyte_umap",
    label_col = "major_celltype",
    target_celltype = "Hepatocyte",
    target_tissue = "Tumor",
    background_color = "#CFCFCF",
    highlight_color = "#FF2D7A",
    contour_color = "#7A7A7A",
    background_size = 0.08,
    highlight_size = 0.16,
    background_alpha = 0.24,
    highlight_alpha = 0.9,
    contour_bins = 13,
    contour_cells = 60000,
    contour_grid = 220,
    seed = 20260601,
    max_cells = 120000,
    width = 6.4,
    height = 6.0,
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
    numeric_args <- c(
      "background_size", "highlight_size", "background_alpha", "highlight_alpha",
      "contour_bins", "contour_cells", "contour_grid", "seed", "max_cells", "width", "height", "dpi"
    )
    if (name %in% numeric_args) {
      value <- as.numeric(value)
    }
    args[[name]] <- value
    i <- i + 2
  }
  args
}

read_table <- function(path, select = NULL) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, select = select)
  } else {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    out <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(select)) {
      out <- out[, intersect(select, names(out)), drop = FALSE]
    }
  }
  first_name <- names(out)[[1]]
  if (first_name == "" || first_name == "V1") {
    names(out)[[1]] <- "cell_id"
  }
  out
}

derive_source_sample <- function(cell_id, dataset, sample_id, study_sample) {
  out <- sample_id
  is_gse149614 <- dataset == "GSE149614"
  prefix <- paste0(study_sample[is_gse149614], "__")
  rest <- cell_id[is_gse149614]
  has_prefix <- startsWith(rest, prefix)
  rest[has_prefix] <- substring(rest[has_prefix], nchar(prefix[has_prefix]) + 1)
  out[is_gse149614] <- sub("_.*$", "", rest)
  out
}

first_nonempty <- function(...) {
  values <- list(...)
  out <- rep("", length(values[[1]]))
  for (value in values) {
    value <- as.character(value)
    fill <- out == "" & !is.na(value) & value != ""
    out[fill] <- value[fill]
  }
  out
}

standardize_tissue <- function(dataset, source_sample, sample_tissue) {
  text <- tolower(as.character(sample_tissue))
  sample <- as.character(source_sample)
  out <- rep("Unknown", length(sample))

  out[grepl("tumou?r", text) | grepl("(^|_)tst$", tolower(sample)) | grepl("^HCC[0-9]+T$", sample)] <- "Tumor"
  out[grepl("normal", text) | grepl("(^|_)nst$", tolower(sample)) | grepl("^HCC[0-9]+N$", sample)] <- "Normal"
  out[grepl("pvtt|portal", text) | grepl("^HCC[0-9]+P$", sample)] <- "PVTT"
  out[grepl("lymph", text) | grepl("^HCC[0-9]+L$", sample)] <- "Lymph"
  out[grepl("cirrh|cst", text) | grepl("(^|_)cst$", tolower(sample))] <- "Cirrhotic"
  out[grepl("single|^sc$", text) | grepl("(^|_)sc$", tolower(sample))] <- "Single_cell_sample"
  out
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

density_contour_grid <- function(df, max_cells, n_grid, seed) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("R package MASS is required to compute density contours.", call. = FALSE)
  }
  density_df <- df[, c("UMAP_1", "UMAP_2"), drop = FALSE]
  density_df <- density_df[stats::complete.cases(density_df), , drop = FALSE]
  if (nrow(density_df) < 10) {
    stop("Too few cells to compute density contours.", call. = FALSE)
  }
  if (is.finite(max_cells) && nrow(density_df) > max_cells) {
    set.seed(as.integer(seed))
    density_df <- density_df[sample(seq_len(nrow(density_df)), size = max_cells, replace = FALSE), , drop = FALSE]
  }
  kde <- MASS::kde2d(density_df$UMAP_1, density_df$UMAP_2, n = as.integer(n_grid))
  contour_df <- expand.grid(UMAP_1 = kde$x, UMAP_2 = kde$y)
  contour_df$density <- as.vector(kde$z)
  contour_df
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(as.integer(args$seed))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(args$plot_data)) {
  message("READ cached plot data ", args$plot_data)
  plot_df <- read_table(args$plot_data)
  if (!"is_target" %in% names(plot_df)) {
    stop("Cached plot data must include is_target.", call. = FALSE)
  }
  plot_df$is_target <- plot_df$is_target %in% c(TRUE, "TRUE", "True", "true", 1, "1")
} else {
  message("BUILD plot data cache ", args$plot_data)
  dir.create(dirname(args$plot_data), recursive = TRUE, showWarnings = FALSE)

message("READ ", args$umap)
umap <- read_table(args$umap)
message("READ ", args$annotations)
annot <- read_table(args$annotations)
message("READ ", args$sample_summary)
sample_summary <- read_table(args$sample_summary)

required_umap <- c("cell_id", "UMAP_1", "UMAP_2")
missing_umap <- setdiff(required_umap, names(umap))
if (length(missing_umap) > 0) {
  stop("UMAP table is missing required columns: ", paste(missing_umap, collapse = ", "), call. = FALSE)
}

required_annot <- c("cell_id", "dataset", "sample_id", "study_sample", args$label_col)
missing_annot <- setdiff(required_annot, names(annot))
if (length(missing_annot) > 0) {
  stop("Annotation table is missing required columns: ", paste(missing_annot, collapse = ", "), call. = FALSE)
}

if (!all(c("dataset", "sample") %in% names(sample_summary))) {
  stop("Sample summary must include dataset and sample columns.", call. = FALSE)
}

sample_tissue_cols <- intersect(
  c("site_values", "Tissue_values", "condition_values", "Disease.status_values", "Type_values"),
  names(sample_summary)
)
if (length(sample_tissue_cols) == 0) {
  stop("Sample summary has no recognized tissue/source columns.", call. = FALSE)
}

sample_summary$sample_tissue_raw <- do.call(first_nonempty, sample_summary[sample_tissue_cols])
sample_summary <- sample_summary[, c("dataset", "sample", "sample_tissue_raw"), drop = FALSE]
names(sample_summary)[names(sample_summary) == "sample"] <- "source_sample"
sample_summary <- sample_summary[!duplicated(sample_summary[, c("dataset", "source_sample")]), , drop = FALSE]

annot$source_sample <- derive_source_sample(
  annot$cell_id,
  annot$dataset,
  annot$sample_id,
  annot$study_sample
)
message("JOIN cell annotations to sample tissue metadata")
if (requireNamespace("data.table", quietly = TRUE)) {
  annot <- data.table::as.data.table(annot)
  sample_summary <- data.table::as.data.table(sample_summary)
  annot <- data.table::merge.data.table(annot, sample_summary, by = c("dataset", "source_sample"), all.x = TRUE, sort = FALSE)
  annot <- as.data.frame(annot)
} else {
  annot <- merge(annot, sample_summary, by = c("dataset", "source_sample"), all.x = TRUE, sort = FALSE)
}
annot$tissue_class <- standardize_tissue(annot$dataset, annot$source_sample, annot$sample_tissue_raw)
annot$is_target <- annot[[args$label_col]] == args$target_celltype & annot$tissue_class == args$target_tissue

message("JOIN UMAP coordinates to target labels")
if (requireNamespace("data.table", quietly = TRUE)) {
  umap_dt <- data.table::as.data.table(umap[, required_umap, drop = FALSE])
  annot_dt <- data.table::as.data.table(annot[, c("cell_id", args$label_col, "dataset", "source_sample", "sample_tissue_raw", "tissue_class", "is_target"), drop = FALSE])
  plot_df <- data.table::merge.data.table(umap_dt, annot_dt, by = "cell_id", all.x = TRUE, sort = FALSE)
  plot_df <- as.data.frame(plot_df)
} else {
  plot_df <- merge(
    umap[, required_umap, drop = FALSE],
    annot[, c("cell_id", args$label_col, "dataset", "source_sample", "sample_tissue_raw", "tissue_class", "is_target"), drop = FALSE],
    by = "cell_id",
    all.x = TRUE,
    sort = FALSE
  )
}
plot_df$is_target[is.na(plot_df$is_target)] <- FALSE
plot_df <- plot_df[stats::complete.cases(plot_df[, c("UMAP_1", "UMAP_2")]), , drop = FALSE]
  con <- gzfile(args$plot_data, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(plot_df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  message("WROTE cached plot data ", args$plot_data)
}

if (is.finite(args$max_cells) && nrow(plot_df) > args$max_cells) {
  keep <- sample(seq_len(nrow(plot_df)), size = args$max_cells, replace = FALSE)
  target_idx <- which(plot_df$is_target)
  keep <- sort(unique(c(keep, target_idx)))
  plot_df <- plot_df[keep, , drop = FALSE]
}

target_df <- plot_df[plot_df$is_target, , drop = FALSE]
if (nrow(target_df) == 0) {
  stop("No target cells found for tissue=", args$target_tissue, " and celltype=", args$target_celltype, call. = FALSE)
}
arrows <- axis_arrows(plot_df)
message("COMPUTE density contours")
contours <- density_contour_grid(plot_df, args$contour_cells, args$contour_grid, args$seed)

message("PLOT")
p <- ggplot() +
  geom_point(
    data = plot_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = args$background_color,
    size = args$background_size,
    alpha = args$background_alpha,
    stroke = 0,
    shape = 16
  ) +
  geom_contour(
    data = contours,
    aes(x = UMAP_1, y = UMAP_2, z = density),
    color = args$contour_color,
    linewidth = 0.22,
    bins = as.integer(args$contour_bins),
    alpha = 0.75
  ) +
  geom_point(
    data = target_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = args$highlight_color,
    size = args$highlight_size,
    alpha = args$highlight_alpha,
    stroke = 0,
    shape = 16
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
  theme(plot.margin = margin(6, 6, 6, 6))

png_path <- file.path(args$output_dir, paste0(args$output_prefix, ".png"))
pdf_path <- file.path(args$output_dir, paste0(args$output_prefix, ".pdf"))
ggsave(png_path, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
ggsave(pdf_path, p, width = args$width, height = args$height, bg = "white", useDingbats = FALSE)

summary_path <- file.path(args$output_dir, paste0(args$output_prefix, "_summary.tsv"))
tissue_counts <- as.data.frame(table(plot_df$tissue_class, plot_df[[args$label_col]]), stringsAsFactors = FALSE)
names(tissue_counts) <- c("tissue_class", "cell_type", "n_cells")
tissue_counts <- tissue_counts[tissue_counts$n_cells > 0, , drop = FALSE]
tissue_counts$is_target_combination <- tissue_counts$tissue_class == args$target_tissue & tissue_counts$cell_type == args$target_celltype
write.table(tissue_counts, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("TARGET cells=", nrow(target_df), " total_cells=", nrow(plot_df))
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
