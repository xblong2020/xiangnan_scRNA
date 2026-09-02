#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(patchwork))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
objects <- c(
  file.path(source_paths$processed, "plot_objects", "figure5a.rds"),
  file.path(paths$processed, "plot_objects", paste0("figure5", letters[2:8], ".rds"))
)
if (any(!file.exists(objects))) stop("Missing panel plot objects: ", paste(basename(objects)[!file.exists(objects)], collapse = ", "))
plots <- lapply(objects, readRDS)
preview <- plots[[1]] / (plots[[2]] | plots[[3]]) / (plots[[4]] | plots[[5]]) / (plots[[6]] | plots[[7]]) / plots[[8]] +
  plot_layout(heights = c(0.78, 1.15, 0.9, 1.05, 0.72)) + plot_annotation(
    title = "Corrected temporal landmarks support overlapping regulatory programmes with partially unresolved strict ordering",
    theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 12, hjust = 0.5))
  )
dir.create(paths$preview, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(paths$preview, "figure5_temporal_positioning_onset_fix_a_to_h_preview.pdf"), preview, device = cairo_pdf,
       width = 13.5, height = 22, units = "in", limitsize = FALSE)
ggsave(file.path(paths$preview, "figure5_temporal_positioning_onset_fix_a_to_h_preview.png"), preview, device = "png", dpi = 300,
       width = 13.5, height = 22, units = "in", bg = "white", limitsize = FALSE)
