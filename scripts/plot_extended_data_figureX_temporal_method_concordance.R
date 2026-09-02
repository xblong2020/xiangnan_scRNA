#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_onset_fix_paths(root, create = FALSE)
dt <- fread(file.path(source_paths$metadata, "figure5f_method_concordance.tsv"))
dt <- dt[!is.na(method) & !is.na(comparison)]
method_order <- unique(dt$method)
dt[, method := factor(method, levels = rev(method_order))]
dt[, comparison := factor(comparison, levels = c("A before B", "B before C", "A before C"))]
dt[, label := ifelse(is.finite(probability), sprintf("%.2f\nn=%d", probability, n_valid),
                     ifelse(grepl("small[- ]n", evidence_weight, ignore.case = TRUE), "small-n\nNot available", "Not available"))]
plot_title <- paste0("Extended Data Figure X\n",
                     "Availability and sensitivity of temporal-order estimates across\n",
                     "trajectory methods and resampling schemes")

p <- ggplot(dt, aes(comparison, method, fill = status)) +
  geom_tile(colour = lancet_palette[8], linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.05, lineheight = 0.9) +
  scale_fill_manual(values = evidence_palette, drop = FALSE) +
  labs(
    title = plot_title,
    subtitle = "All NA, small-n and Not available results are retained; method-specific rows are sensitivity/coverage audits",
    x = NULL, y = NULL, fill = "Evidence"
  ) +
  theme_figure5(7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.text.y = element_text(size = 6.5),
        legend.position = "bottom", plot.title = element_text(size = 9.2, face = "bold"))

out_dir <- file.path(paths$extended_figures, "figureX_temporal_method_concordance")
outputs <- export_figure5_plot(p, file.path(out_dir, "extended_data_figureX_temporal_method_concordance"),
                               9.4, max(6.0, uniqueN(dt$method) * 0.27))
figure5_write_tsv(dt, file.path(paths$extended_metadata, "extended_data_figureX_temporal_method_concordance.tsv"))
figure5_write_json(list(
  figure = "Extended Data Figure X",
  title = "Availability and sensitivity of temporal-order estimates across trajectory methods and resampling schemes",
  retained_statuses = c("NA", "small-n", "Not available"),
  source_namespace = "figure5_temporal_positioning_onset_fix",
  n_rows = nrow(dt),
  outputs = as.list(outputs)
), file.path(paths$extended_metadata, "extended_data_figureX_report.json"))
dir.create(paths$extended_processed, recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$extended_processed, "extended_data_figureX_temporal_method_concordance.rds"))
