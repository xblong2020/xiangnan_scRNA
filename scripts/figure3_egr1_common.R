#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

figure3_project_root <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_arg[grepl("^--file=", file_arg)]
  script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
  root <- normalizePath(
    if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
    mustWork = FALSE
  )
  if (!dir.exists(file.path(root, "scripts"))) root <- normalizePath(getwd(), mustWork = TRUE)
  root
}

figure3_get_arg <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

figure3_state_order <- c(
  "normal_reference", "stressed_injured", "regenerative_progenitor",
  "proliferating_candidate", "malignant_or_malignant_like"
)
figure3_state_labels <- c(
  normal_reference = "Normal/reference",
  stressed_injured = "Stressed/injured",
  regenerative_progenitor = "Regenerative/progenitor",
  proliferating_candidate = "Proliferating candidate",
  malignant_or_malignant_like = "Malignant/malignant-like"
)
figure3_state_palette <- c(
  normal_reference = "#B8B8B8",
  stressed_injured = "#56B4E9",
  regenerative_progenitor = "#009E73",
  proliferating_candidate = "#E69F00",
  malignant_or_malignant_like = "#D55E00"
)
figure3_diverging <- c(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426")

figure3_theme <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
      plot.tag = element_text(size = 10, face = "bold", family = "sans"),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.4),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.34, "cm"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

figure3_save <- function(plot, figure_dir, stem, width, height, tiff = FALSE) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height,
         units = "in", device = grDevices::svg)
  if (tiff) {
    ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
           units = "in", dpi = 600, compression = "lzw")
  }
}

figure3_norm_path <- function(path) enc2utf8(gsub("\\\\", "/", normalizePath(path, mustWork = FALSE)))

figure3_write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_json(value, path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 16)
}
