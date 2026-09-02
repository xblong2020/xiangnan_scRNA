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

control_palette <- c(
  foxm1_cebpb_reference = lancet_palette[5],
  proliferation_control = lancet_palette[6],
  calibration_control = lancet_palette[4],
  random_signature = lancet_palette[8]
)

cohort_palette <- c(
  TCGA_LIHC = lancet_palette[1],
  ICGC_LIRI_JP = lancet_palette[2],
  external_bulk_1 = lancet_palette[3],
  external_bulk_2 = lancet_palette[4],
  external_bulk_3 = lancet_palette[5]
)

evidence_palette <- c(
  Robust = lancet_palette[3],
  Partial = lancet_palette[1],
  Unresolved = lancet_palette[6],
  `Control-dominated` = lancet_palette[2],
  `Not estimable` = "#9E9E9E"
)

figure7_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, colour = "#444444"),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5, colour = "#222222"),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 7),
      axis.line = element_line(linewidth = 0.4, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      plot.margin = margin(6, 9, 6, 6)
    )
}

theme_set(figure7_theme())

figure7_export <- function(plot, stem, width, height, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height,
         units = "in", dpi = dpi, device = "png", bg = "white")
  ggsave(paste0(stem, ".svg"), plot = plot, width = width, height = height,
         units = "in", device = grDevices::svg, bg = "white")
  ggsave(paste0(stem, ".tiff"), plot = plot, width = width, height = height,
         units = "in", dpi = dpi, device = "tiff", compression = "lzw", bg = "white")
  invisible(paste0(stem, c(".pdf", ".png", ".svg", ".tiff")))
}

figure7_panel_label <- function(label) {
  labs(tag = label) + theme(plot.tag = element_text(size = 11, face = "bold"))
}
