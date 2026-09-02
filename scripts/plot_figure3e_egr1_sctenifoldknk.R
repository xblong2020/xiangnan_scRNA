#!/usr/bin/env Rscript

## Figure 3E: FDR-significant genes after EGR1 virtual knockout.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(scales)
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
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3e_egr1")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3e_egr1")),
  mustWork = FALSE
)
subset <- figure3_get_arg("--subset", "stressed_regenerative")
input_path <- normalizePath(
  figure3_get_arg(
    "--input",
    file.path(data_dir, paste0("figure3e_egr1_", subset, "_consensus_perturbation_genes.tsv"))
  ),
  mustWork = TRUE
)
run_report_path <- normalizePath(
  figure3_get_arg(
    "--run-report",
    file.path(data_dir, paste0("figure3e_egr1_", subset, "_run_report.json"))
  ),
  mustWork = TRUE
)
top_n <- as.integer(figure3_get_arg("--top-n", "20"))
fdr_cutoff <- as.numeric(figure3_get_arg("--fdr-cutoff", "0.05"))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("tf", "gene", "distance", "p.adj", "p.value", "Z", "FC", "subset")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing required Figure 3E columns: ", paste(missing, collapse = ", "))
dat$distance <- as.numeric(dat$distance)
dat$p.adj <- as.numeric(dat$p.adj)
dat$p.value <- as.numeric(dat$p.value)
dat$Z <- as.numeric(dat$Z)
dat$FC <- as.numeric(dat$FC)
tf_dat <- dat[
  dat$tf == "EGR1" & dat$gene != "EGR1" & is.finite(dat$distance),
  ,
  drop = FALSE
]
significant <- tf_dat[
  is.finite(tf_dat$p.adj) & tf_dat$p.adj < fdr_cutoff,
  ,
  drop = FALSE
]
significant <- significant[
  order(significant$distance, decreasing = TRUE, na.last = NA),
  ,
  drop = FALSE
]
plot_data <- head(significant, top_n)
plot_data$rank <- seq_len(nrow(plot_data))
plot_data$minus_log10_fdr <- -log10(pmax(plot_data$p.adj, .Machine$double.xmin))
plot_data$significance <- ifelse(
  plot_data$p.adj < 0.001,
  "***",
  ifelse(plot_data$p.adj < 0.01, "**", "*")
)

source_path <- file.path(data_dir, "figure3e_egr1_significant_perturbed_genes.tsv")
output_columns <- c(
  "rank", "tf", "gene", "distance", "Z", "FC", "p.value", "p.adj",
  "minus_log10_fdr", "significance", "subset"
)
if (nrow(plot_data)) {
  write.table(
    plot_data[, output_columns, drop = FALSE],
    source_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else {
  empty <- as.data.frame(setNames(replicate(length(output_columns), logical(0), simplify = FALSE), output_columns))
  write.table(empty, source_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

stem <- file.path(figure_dir, "figure3e_egr1_significant_perturbed_genes")
plot_generated <- nrow(plot_data) > 0
if (plot_generated) {
  plot_data$gene_plot <- factor(plot_data$gene, levels = rev(plot_data$gene))
  lancet <- ggsci::pal_lancet("lanonc")(9)
  gradient_colours <- c(lancet[1], lancet[4], lancet[3], lancet[6], lancet[2])
  x_max <- max(plot_data$distance, na.rm = TRUE)
  title <- if (nrow(significant) >= 20L) {
    "Top 20 significantly perturbed genes"
  } else {
    paste0("Significantly perturbed genes after EGR1 knockout (n = ", nrow(significant), ")")
  }
  p <- ggplot(plot_data, aes(x = distance, y = gene_plot)) +
    geom_col(aes(fill = distance), width = 0.68, alpha = 0.90, colour = NA) +
    geom_point(aes(fill = distance), shape = 21, size = 2.25, stroke = 0.45, colour = "white") +
    geom_text(
      aes(label = significance),
      x = plot_data$distance + x_max * 0.025,
      hjust = 0,
      size = 2.8,
      colour = "grey20"
    ) +
    scale_x_continuous(
      labels = scales::label_scientific(digits = 2),
      expand = expansion(mult = c(0, 0.16))
    ) +
    scale_fill_gradientn(
      colours = gradient_colours,
      labels = scales::label_scientific(digits = 2),
      name = "Distance"
    ) +
    labs(x = "Manifold alignment distance", y = NULL, title = title, tag = "Figure 3E") +
    coord_cartesian(clip = "off") +
    figure3_theme() +
    theme(
      axis.text.y = element_text(size = 7.7, colour = "black"),
      legend.position = "right",
      plot.margin = margin(6, 15, 6, 6)
    )
  height <- max(3.0, min(6.2, 2.1 + 0.155 * nrow(plot_data)))
  figure3_save(p, figure_dir, "figure3e_egr1_significant_perturbed_genes", 5.4, height, tiff = TRUE)
}

run_report <- jsonlite::read_json(run_report_path, simplifyVector = TRUE)
review_risks <- run_report$review_risk_flags
if (!plot_generated) {
  review_risks <- c(
    review_risks,
    list(list(
      flag = "no_fdr_significant_perturbed_genes",
      severity = "main_panel_blocking",
      detail = "No non-EGR1 gene passed the conservative across-seed p.adj < 0.05 rule; the formal Figure 3E plot was suppressed."
    ))
  )
}
report <- list(
  module = "Figure 3E",
  method = "scTenifoldKnk EGR1 virtual knockout with conservative across-seed consensus",
  target_tf = "EGR1",
  target_gene_excluded = TRUE,
  subset = subset,
  input = figure3_norm_path(input_path),
  run_report = figure3_norm_path(run_report_path),
  parameters = run_report$parameters,
  consensus_method = run_report$consensus,
  selection = paste0("gene != EGR1 and p.adj < ", fdr_cutoff, ", then descending manifold-alignment distance"),
  n_tested_excluding_target = nrow(tf_dat),
  n_significant_excluding_target = nrow(significant),
  n_plotted = nrow(plot_data),
  non_significant_fill_used = FALSE,
  formal_plot_generated = plot_generated,
  dynamic_title = if (plot_generated) title else NULL,
  top_genes = as.list(as.character(plot_data$gene)),
  source_table = figure3_norm_path(source_path),
  outputs = if (plot_generated) list(
    pdf = figure3_norm_path(paste0(stem, ".pdf")),
    png = figure3_norm_path(paste0(stem, ".png")),
    svg = figure3_norm_path(paste0(stem, ".svg")),
    tiff = figure3_norm_path(paste0(stem, ".tiff"))
  ) else list(),
  review_risk_flags = review_risks,
  caveat = "Manifold alignment distance is unsigned network displacement. This computational perturbation is not experimental validation and is not interpreted as directional gene regulation."
)
figure3_write_json(report, file.path(data_dir, "figure3e_egr1_report.json"))
message(
  "Figure 3E report written; significant genes=", nrow(significant),
  ", formal plot generated=", plot_generated
)
