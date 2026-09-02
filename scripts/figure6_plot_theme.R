#!/usr/bin/env Rscript

## Shared Figure 6 colour, typography and export contract.

figure6_theme_root <- function() {
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

FIGURE6_PROJECT_ROOT <- figure6_theme_root()
FIGURE6_R_LIBRARY <- file.path(
  FIGURE6_PROJECT_ROOT,
  "data", "processed", "driver", "figure6_directional_network", "r_library"
)
if (dir.exists(FIGURE6_R_LIBRARY)) .libPaths(c(FIGURE6_R_LIBRARY, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(jsonlite)
})

lancet_palette <- ggsci::pal_lancet("lanonc")(9)

axis_palette <- c(
  identity_axis = lancet_palette[1],
  stress_axis = lancet_palette[3],
  sox4_axis = lancet_palette[2]
)

neutral_gray <- lancet_palette[8]
dark_text <- lancet_palette[9]
auxiliary_cyan <- lancet_palette[4]
auxiliary_purple <- lancet_palette[5]
auxiliary_light <- lancet_palette[6]

effect_gradient <- c(
  low = lancet_palette[1],
  mid = "#F7F7F7",
  high = lancet_palette[2]
)

evidence_palette <- c(
  strong = lancet_palette[3],
  moderate = lancet_palette[4],
  weak = lancet_palette[6],
  unresolved = lancet_palette[8],
  opposite = lancet_palette[2],
  unavailable = "#FFFFFF"
)

state_palette <- c(
  normal_reference = lancet_palette[8],
  stressed_injured = lancet_palette[4],
  regenerative_progenitor = lancet_palette[3],
  proliferating_candidate = lancet_palette[5],
  malignant_or_malignant_like = lancet_palette[2]
)

figure6_theme <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(size = 10, face = "plain", hjust = 0.5, colour = dark_text),
      plot.subtitle = element_text(size = 8, colour = dark_text, hjust = 0.5),
      plot.caption = element_text(size = 6.8, colour = dark_text, hjust = 0),
      plot.tag = element_text(size = 10, face = "bold", family = "sans", colour = dark_text),
      axis.title = element_text(size = 8.5, colour = dark_text),
      axis.text = element_text(size = 7.5, colour = dark_text),
      axis.line = element_line(linewidth = 0.4, colour = dark_text),
      legend.title = element_text(size = 7.2, colour = dark_text),
      legend.text = element_text(size = 7, colour = dark_text),
      legend.key.height = grid::unit(0.34, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 6)
    )
}

figure6_save <- function(plot, figure_dir, stem, width, height) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in", device = grDevices::cairo_pdf,
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = 600,
    bg = "white", limitsize = FALSE
  )
  grDevices::svg(
    file.path(figure_dir, paste0(stem, ".svg")),
    width = width, height = height, family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".tiff")), plot,
    width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white", limitsize = FALSE
  )
  invisible(file.path(figure_dir, paste0(stem, c(".pdf", ".png", ".svg", ".tiff"))))
}

figure6_palette_contract <- function() {
  list(
    r_version = R.version.string,
    ggsci_version = as.character(utils::packageVersion("ggsci")),
    lancet_palette = as.list(unname(lancet_palette)),
    axis_palette = as.list(axis_palette),
    effect_gradient = as.list(effect_gradient),
    evidence_palette = as.list(evidence_palette),
    midpoint = 0,
    source = "ggsci::pal_lancet('lanonc')(9)"
  )
}

