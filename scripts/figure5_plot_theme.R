#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(scales)
})

lancet_palette <- ggsci::pal_lancet("lanonc")(9)
axis_palette <- c(
  identity_loss = lancet_palette[1],
  stress_transition = lancet_palette[3],
  sox4_stabilization = lancet_palette[2]
)
evidence_palette <- c(
  "Supported" = lancet_palette[3],
  "Partial" = lancet_palette[6],
  "Unstable" = lancet_palette[5],
  "Not resolved" = lancet_palette[8],
  "Opposite" = lancet_palette[2],
  "Not available" = "#FFFFFF"
)
figure5_diverging_palette <- c(
  lancet_palette[1], lancet_palette[4], "#F7F7F7", lancet_palette[6], lancet_palette[2]
)

theme_figure5 <- function(base_size = 7.5) {
  theme_classic(base_family = "sans", base_size = base_size) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, colour = lancet_palette[9]),
      axis.title = element_text(size = 8.5, colour = lancet_palette[9]),
      axis.text = element_text(size = 7.5, colour = lancet_palette[9]),
      axis.line = element_line(linewidth = 0.4, colour = lancet_palette[9]),
      axis.ticks = element_line(linewidth = 0.4, colour = lancet_palette[9]),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = grid::unit(2.5, "mm"),
      plot.margin = margin(5, 7, 5, 7)
    )
}

export_figure5_plot <- function(plot, stem, width, height, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  outputs <- c(
    pdf = paste0(stem, ".pdf"),
    png = paste0(stem, ".png"),
    svg = paste0(stem, ".svg"),
    tiff = paste0(stem, ".tiff")
  )
  ggplot2::ggsave(outputs[["pdf"]], plot = plot, device = grDevices::cairo_pdf,
                  width = width, height = height, units = "in", limitsize = FALSE)
  ggplot2::ggsave(outputs[["png"]], plot = plot, device = "png", dpi = dpi,
                  width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
  grDevices::tiff(outputs[["tiff"]], width = width, height = height, units = "in",
                  res = dpi, compression = "lzw", bg = "white", family = "sans")
  print(plot)
  grDevices::dev.off()
  grDevices::svg(outputs[["svg"]], width = width, height = height, family = "sans", bg = "white")
  print(plot)
  grDevices::dev.off()
  outputs
}
