#!/usr/bin/env Rscript

## Review-only Figure 3E sensitivity summary across state definitions.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "scripts", "figure3_egr1_common.R"))

data_dir <- normalizePath(
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3e_egr1_sensitivity")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3e_egr1_sensitivity")),
  mustWork = FALSE
)
summary <- read.delim(
  file.path(data_dir, "figure3e_egr1_sensitivity_summary.tsv"),
  stringsAsFactors = FALSE
)
pairwise <- read.delim(
  file.path(data_dir, "figure3e_egr1_sensitivity_pairwise.tsv"),
  stringsAsFactors = FALSE
)
subset_labels <- c(
  stressed_injured = "Stressed",
  stressed_regenerative = "Stressed + regenerative",
  intermediate_pseudotime = "Intermediate pseudotime",
  malignant_like = "Malignant-like"
)
summary$subset_label <- factor(
  unname(subset_labels[summary$subset]),
  levels = unname(subset_labels[summary$subset])
)
summary$analysis_role <- ifelse(summary$selected_main, "Main", "Sensitivity")

p_count <- ggplot(
  summary,
  aes(x = subset_label, y = n_significant_perturbed_genes, fill = analysis_role)
) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n_significant_perturbed_genes), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c(Main = "#D55E00", Sensitivity = "#56B4E9"), name = NULL) +
  labs(
    x = NULL,
    y = "FDR-significant perturbed genes",
    title = "EGR1 network sensitivity"
  ) +
  figure3_theme() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  )

pairwise$left_label <- unname(subset_labels[pairwise$subset_left])
pairwise$right_label <- unname(subset_labels[pairwise$subset_right])
p_jaccard <- ggplot(
  pairwise,
  aes(x = left_label, y = right_label, fill = significant_gene_jaccard)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = ifelse(is.finite(significant_gene_jaccard), sprintf("%.2f", significant_gene_jaccard), "NA")), size = 2.6) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#009E73",
    limits = c(0, 1),
    na.value = "grey90",
    name = "Jaccard"
  ) +
  labs(x = NULL, y = NULL, title = "Significant-gene overlap") +
  figure3_theme() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

combined <- p_count + p_jaccard + plot_layout(widths = c(1, 1.15))
figure3_save(
  combined,
  figure_dir,
  "figure3e_egr1_state_sensitivity_summary",
  10.5,
  4.4,
  tiff = TRUE
)
report <- list(
  module = "Figure 3E sensitivity plotting",
  target_tf = "EGR1",
  main_subset = "stressed_regenerative",
  plotting_language = "R",
  r_version = R.version.string,
  outputs = list(
    pdf = figure3_norm_path(file.path(figure_dir, "figure3e_egr1_state_sensitivity_summary.pdf")),
    png = figure3_norm_path(file.path(figure_dir, "figure3e_egr1_state_sensitivity_summary.png")),
    svg = figure3_norm_path(file.path(figure_dir, "figure3e_egr1_state_sensitivity_summary.svg")),
    tiff = figure3_norm_path(file.path(figure_dir, "figure3e_egr1_state_sensitivity_summary.tiff"))
  ),
  caveat = "Review-only sensitivity plot; non-main networks use lower replication and are not substitutes for the formal main analysis."
)
figure3_write_json(report, file.path(data_dir, "figure3e_egr1_sensitivity_r_plot_report.json"))

