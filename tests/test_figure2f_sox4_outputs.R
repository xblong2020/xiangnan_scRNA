#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- args0[grepl("^--file=", args0)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
root <- normalizePath(if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(), mustWork = TRUE)

metadata_dir <- file.path(root, "metadata/driver/figure2f_sox4")
figure_dir <- file.path(root, "figures/driver/figure2f_sox4")
plot_path <- file.path(metadata_dir, "figure2f_sox4_plot_data.tsv")
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 0)
dat <- read.delim(plot_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(dat) >= 1L, nrow(dat) <= 10L)
stopifnot(all(diff(dat$minus_log10_pvalue) <= 0))
stopifnot(all(dat$overlap_count > 0), all(dat$pvalue > 0), all(dat$pvalue <= 1))
stopifnot(all(dat$tf == "SOX4"), all(dat$subset == "malignant_like"))

required <- c(
  file.path(metadata_dir, "figure2f_sox4_enrichment_all.tsv"),
  file.path(metadata_dir, "figure2f_sox4_report.json"),
  paste0(file.path(figure_dir, "figure2f_sox4_pathway_enrichment"), c(".png", ".pdf", ".svg", ".tiff"))
)
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0))
report <- jsonlite::fromJSON(file.path(metadata_dir, "figure2f_sox4_report.json"))
stopifnot(report$n_significant_perturbed_genes == 71L)
stopifnot(report$n_plotted == nrow(dat))
message("Figure 2F output verification passed: ", nrow(dat), " pathways")
