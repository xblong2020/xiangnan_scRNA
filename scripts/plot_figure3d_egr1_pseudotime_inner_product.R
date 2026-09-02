#!/usr/bin/env Rscript

## Figure 3D: EGR1 perturbation score along strict-main pseudotime.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(scales)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "scripts", "figure3_egr1_common.R"))

c_data_dir <- normalizePath(
  figure3_get_arg("--figure3c-data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3c_egr1")),
  mustWork = FALSE
)
data_dir <- normalizePath(
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3d_egr1")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3d_egr1")),
  mustWork = FALSE
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

umap <- read.delim(
  gzfile(file.path(data_dir, "figure3d_egr1_pseudotime_inner_product_umap.tsv.gz")),
  stringsAsFactors = FALSE
)
tsne <- read.delim(
  gzfile(file.path(data_dir, "figure3d_egr1_pseudotime_inner_product_tsne.tsv.gz")),
  stringsAsFactors = FALSE
)
bins <- read.delim(
  file.path(data_dir, "figure3d_egr1_pseudotime_bin_summary.tsv"),
  stringsAsFactors = FALSE
)
statistics_report <- jsonlite::read_json(
  file.path(data_dir, "figure3d_egr1_report.json"),
  simplifyVector = TRUE
)
score_report <- jsonlite::read_json(
  file.path(c_data_dir, "figure3c_egr1_inner_product_report.json"),
  simplifyVector = TRUE
)
shared_limit <- as.numeric(score_report$colour_scale$three_axis_shared_symmetric_limit)
panel_limit <- as.numeric(score_report$colour_scale$panel_specific_symmetric_limit)

plot_relation <- function(df, bin_data, title, score_limit, tag = NULL) {
  sub <- df[df$valid %in% TRUE, , drop = FALSE]
  ggplot(sub, aes(x = pseudotime_grid, y = inner_product_score_grid)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      span = 0.75,
      linewidth = 0.45,
      colour = "grey35",
      fill = "grey82",
      alpha = 0.45
    ) +
    geom_point(aes(colour = inner_product_score_grid), size = 1.45, alpha = 0.78) +
    geom_line(
      data = bin_data[order(bin_data$pseudotime_center), ],
      aes(x = pseudotime_center, y = score_median),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.55
    ) +
    geom_point(
      data = bin_data,
      aes(x = pseudotime_center, y = score_median),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      stroke = 0.45,
      fill = "white",
      colour = "black"
    ) +
    scale_colour_gradient2(
      low = figure3_diverging[["low"]],
      mid = figure3_diverging[["mid"]],
      high = figure3_diverging[["high"]],
      midpoint = 0,
      limits = c(-score_limit, score_limit),
      oob = scales::squish,
      name = "Perturbation score (PS)"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_cartesian(ylim = c(-score_limit, score_limit)) +
    labs(
      x = "driver_main_strict__pseudotime_rank",
      y = "EGR1 perturbation vector · baseline developmental vector",
      title = title,
      tag = tag
    ) +
    figure3_theme() +
    theme(legend.position = "right")
}

umap_bins <- bins[bins$space == "CellOracle UMAP grid", , drop = FALSE]
tsne_bins <- bins[bins$space == "Expression t-SNE grid", , drop = FALSE]

save_set <- function(limit, suffix, official = FALSE) {
  p_umap <- plot_relation(
    umap,
    umap_bins,
    "EGR1 knockout",
    limit,
    if (official) "Figure 3D" else NULL
  )
  p_tsne <- plot_relation(tsne, tsne_bins, "EGR1 knockout (t-SNE sensitivity)", limit)
  stem_suffix <- if (nzchar(suffix)) paste0("_", suffix) else ""
  figure3_save(
    p_umap,
    figure_dir,
    paste0("figure3d_egr1_pseudotime_inner_product_umap", stem_suffix),
    5.7,
    4.3
  )
  figure3_save(
    p_tsne,
    figure_dir,
    paste0("figure3d_egr1_pseudotime_inner_product_tsne", stem_suffix),
    5.7,
    4.3
  )
  combined <- p_umap + p_tsne + plot_layout(guides = "collect") & theme(legend.position = "right")
  figure3_save(
    combined,
    figure_dir,
    paste0("figure3d_egr1_pseudotime_inner_product_umap_tsne", stem_suffix),
    10.8,
    4.3
  )
}

# Main files use the same three-axis shared symmetric limit as Figure 3C.
save_set(shared_limit, "", official = TRUE)
save_set(panel_limit, "panel_specific")

report <- list(
  module = "Figure 3D plotting",
  target_tf = "EGR1",
  plotting_language = "R",
  r_version = R.version.string,
  score_definition = "EGR1 perturbation vector dot baseline developmental vector",
  pseudotime = "driver_main_strict__pseudotime_rank",
  n_bins = 10,
  fixed_stages = list(
    early = "0.00-0.33",
    intermediate = "0.33-0.67",
    late = "0.67-1.00"
  ),
  loess_sensitivity = list(method = "loess", span = 0.75, confidence_interval = 0.95),
  observed_umap_pattern = statistics_report$observed_umap_pattern,
  spearman_bootstrap = statistics_report$spearman_bootstrap,
  colour_scale = list(
    low = unname(figure3_diverging[["low"]]),
    mid = unname(figure3_diverging[["mid"]]),
    high = unname(figure3_diverging[["high"]]),
    midpoint = 0,
    main_shared_limit = shared_limit,
    panel_specific_limit = panel_limit
  ),
  outputs = list(
    umap_pdf = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap.pdf")),
    umap_png = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap.png")),
    umap_svg = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap.svg")),
    tsne_pdf = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_tsne.pdf")),
    tsne_png = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_tsne.png")),
    tsne_svg = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_tsne.svg")),
    combined_pdf = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap_tsne.pdf")),
    combined_png = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap_tsne.png")),
    combined_svg = figure3_norm_path(file.path(figure_dir, "figure3d_egr1_pseudotime_inner_product_umap_tsne.svg"))
  ),
  caveat = "Fixed stages are primary; LOESS, data-driven change points, and t-SNE are sensitivity analyses. The observed sign and peak stage were not constrained."
)
figure3_write_json(report, file.path(data_dir, "figure3d_egr1_r_plot_report.json"))
message("Figure 3D written to: ", figure_dir)

