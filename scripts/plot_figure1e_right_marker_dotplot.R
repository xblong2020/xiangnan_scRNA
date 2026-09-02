suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1e_right_marker_dotplot.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_args <- function(argv) {
  defaults <- list(
    marker_scores = file.path(ROOT, "metadata/hepatocyte/hepatocyte_state_marker_scores.tsv"),
    summary = file.path(ROOT, "metadata/hepatocyte/hepatocyte_subcluster_summary.tsv"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1E_right_marker_dotplot",
    width = 10.5,
    height = 4.8,
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
  malignant_hepatocyte_candidate_needs_cnv = "Malignant candidate (needs CNV)",
  ambiguous_epithelial_or_mixed = "Ambiguous epithelial / mixed"
)

gene_panel <- list(
  Mature_Hepatocyte = c("ALB", "APOA1", "TTR", "HNF4A"),
  Stressed_Injured = c("HSPA1A", "DNAJB1", "FOS", "ATF3"),
  Regenerative_Progenitor = c("EPCAM", "KRT19", "SOX9", "TACSTD2"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "PCNA"),
  HCC_Malignant_Associated = c("GPC3", "AFP", "SPP1", "MDK")
)

read_table <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)

marker_df <- read_table(params$marker_scores)
summary_df <- read_table(params$summary)[, c("leiden_hep", "hepatocyte_state_label", "n_cells"), drop = FALSE]
summary_df$leiden_hep <- as.character(summary_df$leiden_hep)
marker_df$leiden_hep <- as.character(marker_df$leiden_hep)

marker_df <- merge(marker_df, summary_df, by = "leiden_hep", all.x = TRUE, suffixes = c("", "_cluster"))
panel_genes <- unlist(gene_panel, use.names = FALSE)
marker_df <- marker_df[marker_df$gene %in% panel_genes, , drop = FALSE]
gene_to_group <- unlist(
  lapply(
    names(gene_panel),
    function(group_name) stats::setNames(rep(group_name, length(gene_panel[[group_name]])), gene_panel[[group_name]])
  )
)
marker_df$gene_group <- unname(gene_to_group[marker_df$gene])
marker_df <- marker_df[!duplicated(marker_df[, c("leiden_hep", "gene")]), , drop = FALSE]

agg <- aggregate(
  cbind(raw_sum, detected_cells, n_cells_gene_present) ~ hepatocyte_state_label + gene + gene_group,
  data = marker_df,
  FUN = sum
)

weighted_log <- stats::aggregate(
  x = marker_df$log1p_cpm * marker_df$n_cells_gene_present,
  by = list(
    hepatocyte_state_label = marker_df$hepatocyte_state_label,
    gene = marker_df$gene,
    gene_group = marker_df$gene_group
  ),
  FUN = sum
)
names(weighted_log)[names(weighted_log) == "x"] <- "weighted_log1p"
agg <- merge(agg, weighted_log, by = c("hepatocyte_state_label", "gene", "gene_group"), all.x = TRUE)
agg$pct_expr <- agg$detected_cells / pmax(agg$n_cells_gene_present, 1)
agg$mean_log1p_cpm <- agg$weighted_log1p / pmax(agg$n_cells_gene_present, 1)

agg$hepatocyte_state_label <- factor(agg$hepatocyte_state_label, levels = state_order)
gene_order <- unlist(gene_panel, use.names = FALSE)
agg$gene <- factor(agg$gene, levels = gene_order)

agg$scaled_log1p <- 0
for (gene in levels(agg$gene)) {
  idx <- agg$gene == gene
  vals <- agg$mean_log1p_cpm[idx]
  if (all(is.na(vals))) {
    next
  }
  rng <- range(vals, finite = TRUE)
  if (diff(rng) <= 0) {
    agg$scaled_log1p[idx] <- 0.5
  } else {
    agg$scaled_log1p[idx] <- (vals - rng[1]) / diff(rng)
  }
}

group_df <- data.frame(
  gene = factor(gene_order, levels = gene_order),
  gene_group = rep(names(gene_panel), times = lengths(gene_panel)),
  stringsAsFactors = FALSE
)
group_centers <- aggregate(as.numeric(gene) ~ gene_group, data = group_df, FUN = mean)
names(group_centers)[2] <- "xcenter"

separator_x <- cumsum(lengths(gene_panel))[-length(gene_panel)] + 0.5

expr_palette <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
group_palette <- stats::setNames(expr_palette, names(gene_panel))

p <- ggplot(agg, aes(x = gene, y = hepatocyte_state_label)) +
  geom_point(aes(size = pct_expr, fill = scaled_log1p), shape = 21, color = "#5A5A5A", stroke = 0.18) +
  scale_fill_gradientn(colours = expr_palette, limits = c(0, 1), name = "Scaled\nexpression") +
  scale_size(range = c(0.8, 7.2), limits = c(0, 1), name = "Fraction\nexpressed") +
  geom_vline(xintercept = separator_x, color = "#D0D0D0", linewidth = 0.4) +
  scale_y_discrete(labels = state_labels) +
  coord_cartesian(ylim = c(0.6, length(state_order) + 0.7), clip = "off") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEFEF", linewidth = 0.35),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(22, 10, 10, 10),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p <- p + geom_text(
  data = group_centers,
  aes(x = xcenter, y = Inf, label = gsub("_", "\n", gene_group, fixed = TRUE)),
  inherit.aes = FALSE,
  vjust = 1.15,
  size = 3.2,
  fontface = "bold",
  color = unname(group_palette[group_centers$gene_group])
)

png_path <- file.path(params$output_dir, paste0(params$output_prefix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(params$output_prefix, ".pdf"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
