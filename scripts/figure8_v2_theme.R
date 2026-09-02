#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

v1_r_library <- figure8_v2_existing_r_library()
.libPaths(c(v1_r_library, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

figure8_v2_lancet <- ggsci::pal_lancet("lanonc")(9)
figure8_v2_axis_palette <- c(
  axis_A_identity = "#00468BFF",
  axis_B_stress = "#42B540FF",
  axis_C_sox4 = "#ED0000FF",
  malignant_fate = figure8_v2_lancet[[5]]
)
figure8_v2_method_palette <- c(
  DrugReflector = figure8_v2_lancet[[5]],
  L1000FWD = figure8_v2_lancet[[4]],
  CLUE = figure8_v2_lancet[[1]],
  PRISM = figure8_v2_lancet[[3]],
  unavailable = "#D9D9D9"
)
figure8_v2_tier_palette <- c(
  tier_A = figure8_v2_lancet[[3]],
  tier_B = figure8_v2_lancet[[4]],
  tier_C = figure8_v2_lancet[[6]],
  discordant = figure8_v2_lancet[[2]],
  unresolved = figure8_v2_lancet[[8]]
)
figure8_v2_diverging <- c(low = figure8_v2_lancet[[2]], mid = "#F7F7F7", high = figure8_v2_lancet[[1]])

figure8_v2_theme <- function(base_size = 6.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.caption = element_text(size = base_size - 1, colour = "#555555", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      legend.key = element_blank()
    )
}

figure8_v2_save_plot <- function(plot, stem, width_mm = 183, height_mm = 230, dpi = 600) {
  figure8_v2_init_dirs()
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  pdf_path <- file.path(FIGURE8_V2_FIGURES, paste0(stem, ".pdf"))
  svg_path <- file.path(FIGURE8_V2_FIGURES, paste0(stem, ".svg"))
  png_path <- file.path(FIGURE8_V2_FIGURES, paste0(stem, ".png"))
  tiff_path <- file.path(FIGURE8_V2_FIGURES, paste0(stem, ".tiff"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "sans", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  grDevices::svg(svg_path, width = width, height = height, family = "sans", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi, type = "cairo", bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::tiff(tiff_path, width = width, height = height, units = "in", res = dpi, compression = "lzw", type = "cairo", bg = "white")
  print(plot)
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, svg = svg_path, png = png_path, tiff = tiff_path))
}
