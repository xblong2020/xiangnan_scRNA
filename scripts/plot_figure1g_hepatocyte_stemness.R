suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1g_hepatocyte_stemness.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
source(file.path(ROOT, "scripts", "figure1g_stemness_helpers.R"))

parse_args <- function(argv) {
  defaults <- list(
    scores = file.path(ROOT, "metadata/figure1c/figure1c_cytotrace2_scores_by_cell.hepatocyte.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1G_hepatocyte_lineage_cytotrace2_stemness",
    max_points_per_state = 450,
    seed = 20260710,
    width = 7.2,
    height = 5.0,
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
    if (name %in% c("max_points_per_state", "seed", "width", "height", "dpi")) value <- as.numeric(value)
    args[[name]] <- value
    i <- i + 2L
  }
  args
}

read_table <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

sample_points_by_state <- function(df, max_points_per_state, seed) {
  set.seed(as.integer(seed))
  by_state <- split(df, df$hepatocyte_state_label, drop = TRUE)
  sampled <- lapply(by_state, function(state_df) {
    keep_n <- min(nrow(state_df), as.integer(max_points_per_state))
    state_df[sample(seq_len(nrow(state_df)), size = keep_n, replace = FALSE), , drop = FALSE]
  })
  do.call(rbind, sampled)
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
  normal_hepatocyte_like = "Normal\nhepatocyte-like",
  stressed_injured_hepatocyte = "Stressed /\ninjured",
  regenerative_progenitor_like_hepatocyte = "Regenerative /\nprogenitor-like",
  proliferating_hepatocyte_candidate = "Proliferating\ncandidate",
  malignant_hepatocyte_candidate_needs_cnv = "Malignant candidate\n(needs CNV)",
  ambiguous_epithelial_or_mixed = "Ambiguous epithelial /\nmixed"
)

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_data <- prepare_figure1g_data(read_table(params$scores), state_order)
if (nrow(plot_data) == 0L) stop("No finite CytoTRACE2 scores matched the Figure 1E hepatocyte states.", call. = FALSE)

summary_df <- summarize_figure1g_stemness(plot_data, state_order)
point_data <- sample_points_by_state(plot_data, params$max_points_per_state, params$seed)
state_palette <- setNames(ggsci::pal_lancet()(length(state_order)), state_order)

p <- ggplot(plot_data, aes(x = hepatocyte_state_label, y = CytoTRACE2_Score, fill = hepatocyte_state_label)) +
  geom_violin(trim = FALSE, scale = "width", color = "#333333", linewidth = 0.28, alpha = 0.82) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "#333333", linewidth = 0.32) +
  geom_point(
    data = point_data,
    aes(color = hepatocyte_state_label),
    position = position_jitter(width = 0.13, height = 0),
    size = 0.48,
    alpha = 0.34,
    inherit.aes = TRUE
  ) +
  scale_fill_manual(values = state_palette, guide = "none") +
  scale_color_manual(values = state_palette, guide = "none") +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = "CytoTRACE2 stemness score",
    title = "Figure 1G. Stemness scores across hepatocyte lineage states",
    subtitle = "Violin width shows cell density; boxplots show median and interquartile range"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.line.x = element_line(linewidth = 0.35),
    axis.line.y = element_line(linewidth = 0.35),
    axis.text.x = element_text(size = 7.2, color = "black", margin = margin(t = 5)),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8, hjust = 0, margin = margin(b = 6)),
    plot.margin = margin(8, 8, 6, 8)
  )

png_path <- file.path(params$output_dir, paste0(params$output_prefix, ".png"))
pdf_path <- file.path(params$output_dir, paste0(params$output_prefix, ".pdf"))
summary_path <- file.path(params$output_dir, paste0(params$output_prefix, "_summary.tsv"))
ggsave(png_path, p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(pdf_path, p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("WROTE ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
