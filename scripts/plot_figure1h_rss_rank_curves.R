suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/plot_figure1h_rss_rank_curves.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "metadata"))) {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
source(file.path(ROOT, "scripts", "figure1h_rank_curve_helpers.R"))

parse_args <- function(argv) {
  defaults <- list(
    input = file.path(ROOT, "figures/figure1/figure1H_all_hepatocyte_states_rss_source_all_states.tsv.gz"),
    output_dir = file.path(ROOT, "figures/figure1"),
    output_prefix = "figure1H_all_hepatocyte_states_rss_rank_curves",
    top_n = 5,
    width = 17.0,
    height = 4.6,
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
    if (name %in% c("top_n", "width", "height", "dpi")) value <- as.numeric(value)
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
  ambiguous_epithelial_or_mixed = "Ambiguous epithelial /\nmixed (n = 8)"
)

params <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
ranked <- prepare_rss_rank_curves(read_table(params$input), state_order, params$top_n)
ranked$state_display <- factor(state_labels[as.character(ranked$state)], levels = unname(state_labels[state_order]))
state_palette <- setNames(lancet_state_palette(length(state_order)), unname(state_labels[state_order]))
top_df <- build_top_label_positions(ranked)
top_df$state_display <- factor(state_labels[as.character(top_df$state)], levels = unname(state_labels[state_order]))

p <- ggplot(ranked, aes(x = rank, y = rss, group = state_display, color = state_display)) +
  geom_line(linewidth = 0.85, lineend = "round") +
  geom_point(size = 0.72, stroke = 0, alpha = 0.95) +
  geom_segment(
    data = top_df,
    aes(x = point_x, y = point_y, xend = label_x, yend = label_y),
    inherit.aes = FALSE,
    color = "#ED0000FF",
    linewidth = 0.28,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_text(
    data = top_df,
    aes(x = label_x, y = label_y, label = top_label),
    inherit.aes = FALSE,
    color = "#ED0000FF",
    angle = 45,
    hjust = 0,
    vjust = 0,
    size = 2.15,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~state_display, nrow = 1, scales = "free_y") +
  scale_color_manual(values = state_palette, guide = "none", drop = FALSE) +
  scale_x_continuous(
    breaks = c(1, 10, 20, 30),
    limits = c(1, max(ranked$rank) + 1.5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.38))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Figure 1H. Regulon specificity score ranking",
    subtitle = "Regulons are ordered by decreasing RSS within each hepatocyte lineage state; Top 5 are labelled",
    x = "Regulon rank",
    y = "RSS"
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8, color = "#555555", margin = margin(b = 7)),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
    strip.text = element_text(size = 7.5, face = "bold", margin = margin(4, 2, 4, 2)),
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "#D6D6D6", linewidth = 0.35),
    panel.grid.minor = element_line(color = "#ECECEC", linewidth = 0.25),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.spacing.x = grid::unit(0.22, "cm"),
    plot.margin = margin(8, 8, 7, 8)
  )

base <- file.path(params$output_dir, params$output_prefix)
ggsave(paste0(base, ".png"), p, width = params$width, height = params$height, dpi = params$dpi, bg = "white")
ggsave(paste0(base, ".pdf"), p, width = params$width, height = params$height, bg = "white", useDingbats = FALSE)
write.table(top_df, paste0(base, "_top5_regulons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("WROTE ", normalizePath(paste0(base, ".png"), winslash = "/", mustWork = FALSE))
message("WROTE ", normalizePath(paste0(base, ".pdf"), winslash = "/", mustWork = FALSE))
