#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(patchwork))

paths <- figure5_six_panel_paths(root)
obj_dir <- file.path(paths$processed, "plot_objects")
required <- file.path(obj_dir, paste0("figure5_panel_", LETTERS[1:6], ".rds"))
if (!all(file.exists(required))) stop("Run panel A-F scripts before montage", call. = FALSE)
plots <- lapply(required, readRDS)
names(plots) <- LETTERS[1:6]

caption <- paste(
  "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events.",
  "Because eligible patient-token data did not span the complete pseudotemporal range, the primary bootstrap used the prespecified sample-token coverage fallback.",
  sep = "\n"
)

montage <- plots$A /
  (plots$B | plots$C) /
  (plots$D | plots$E) /
  plots$F +
  plot_layout(heights = c(1.12, 1.85, 1.15, 1.22), guides = "keep") +
  plot_annotation(
    title = "Corrected temporal positioning reveals overlapping regulatory programmes with later SOX4 prominence.",
    subtitle = "Figure 5. Coverage-qualified frozen programmes, corrected bootstrap landmarks and conservative activity bands",
    caption = caption,
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
                  plot.subtitle = element_text(size = 8.5, colour = lancet_palette[8], hjust = 0),
                  plot.caption = element_text(size = 7, colour = lancet_palette[9], hjust = 0, lineheight = 0.95),
                  plot.margin = margin(5, 8, 5, 8))
  )

outputs <- export_figure5_plot(montage, file.path(paths$preview, "figure5_six_panel_main"), 13.5, 18.2)
saveRDS(montage, file.path(paths$processed, "figure5_six_panel_montage.rds"))
figure5_write_json(list(
  figure = "Figure 5",
  title = "Corrected temporal positioning reveals overlapping regulatory programmes with later SOX4 prominence.",
  panels = c(
    "5A Corrected temporal-positioning analysis workflow",
    "5B Coverage-corrected regulatory programmes along unified pseudotime",
    "5C Temporal organization of frozen TF, regulon and target-gene programmes",
    "5D Corrected bootstrap temporal landmarks",
    "5E Tie-aware temporal precedence across complementary landmarks",
    "5F Conservative overlapping regulatory-activity model"
  ),
  layout = "row1: A full width; row2: B + C; row3: D + E; row4: F full width",
  caption = caption,
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_montage_report.json"))

