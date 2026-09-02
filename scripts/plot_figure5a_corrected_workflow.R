#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(ggplot2))

paths <- figure5_six_panel_paths(root)

workflow <- data.table::data.table(
  step = 1:9,
  title = c(
    "Frozen hepatocyte atlas\nand metadata",
    "Frozen three-axis\nprogrammes:\n• HNF4A/PPARA\n  identity loss\n• AP-1/CEBPB/EGR1\n  stress transition\n• SOX4 malignant-state\n  stabilization",
    "Oriented consensus\npseudotime:\n0 = normal/reference\n1 = CNV-supported\nmalignant/malignant-like",
    "CNV-supported\nmalignant endpoint\nand CellRank\nmalignant-fate\nprobability",
    "Coverage-qualified\npatient/sample-token\npseudobulk GAM",
    "Corrected temporal\nlandmarks:\nonset • t10 • t50\n• maximum slope",
    "Dataset-stratified\nbootstrap",
    "Tie-aware\nprecedence\nprobability",
    "Conservative\noverlapping\nregulatory-activity\nmodel"
  ),
  fill = c(lancet_palette[4], lancet_palette[1], lancet_palette[3], lancet_palette[6],
           lancet_palette[2], lancet_palette[5], lancet_palette[7], lancet_palette[8], lancet_palette[9])
)
workflow[, `:=`(xmin = step - 0.46, xmax = step + 0.46, ymin = 1.98, ymax = 3.62)]

state_labels <- data.table::data.table(
  state = c("Normal/reference", "Stressed/injured", "Regenerative/progenitor",
            "Proliferating candidate", "Malignant/malignant-like"),
  x = seq(1, 9, length.out = 5)
)
workflow_formal <- c(
  "Frozen hepatocyte atlas and metadata",
  "Frozen three-axis programmes: HNF4A/PPARA identity loss; AP-1/CEBPB/EGR1 stress transition; SOX4 malignant-state stabilization",
  "Oriented consensus pseudotime: 0 = normal/reference; 1 = CNV-supported malignant/malignant-like",
  "CNV-supported malignant endpoint and CellRank malignant-fate probability",
  "Coverage-qualified patient/sample-token pseudobulk GAM",
  "Corrected temporal landmarks: onset; t10; t50; maximum slope",
  "Dataset-stratified bootstrap",
  "Tie-aware precedence probability",
  "Conservative overlapping regulatory-activity model"
)

p <- ggplot() +
  geom_rect(data = workflow,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = "white", linewidth = 0.8) +
  geom_text(data = workflow, aes(x = step, y = (ymin + ymax) / 2, label = title),
            colour = "white", fontface = "bold", size = 1.82, lineheight = 0.84) +
  geom_segment(data = workflow[-nrow(workflow)],
               aes(x = xmax + 0.02, xend = xmin + 0.98, y = 2.80, yend = 2.80),
               arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed"),
               linewidth = 0.55, colour = lancet_palette[9]) +
  annotate("text", x = 5, y = 1.62, label = "Descriptive state anchors along oriented consensus pseudotime",
           size = 2.65, fontface = "bold", colour = lancet_palette[9]) +
  geom_segment(data = state_labels, aes(x = x, xend = x, y = 1.42, yend = 1.26),
               linewidth = 0.45, linetype = "dotted", colour = lancet_palette[8]) +
  geom_point(data = state_labels, aes(x = x, y = 1.20), size = 2.2,
             colour = lancet_palette[3], fill = "white", shape = 21, stroke = 0.7) +
  geom_text(data = state_labels, aes(x = x, y = 0.94, label = state),
            size = 2.15, lineheight = 0.9, colour = lancet_palette[9]) +
  annotate("text", x = 5, y = 0.42,
           label = "State labels provide descriptive anchors and do not define discrete temporal transitions.",
           size = 2.35, fontface = "italic", colour = lancet_palette[9]) +
  annotate("text", x = 5, y = 0.15,
           label = "R / ggplot2 / ggsci::pal_lancet(\"lanonc\") • coverage-qualified, frozen-programme workflow",
           size = 2.05, colour = lancet_palette[8]) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.30, 9.70), breaks = workflow$step,
                     labels = paste0("", workflow$step), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3.88), breaks = NULL, expand = c(0, 0)) +
  labs(title = "A  Corrected temporal-positioning analysis workflow", x = "Analysis sequence", y = NULL,
       subtitle = "Frozen programmes, oriented consensus pseudotime, coverage-qualified GAM and conservative temporal ordering") +
  theme_void(base_family = "sans", base_size = 7.5) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0, colour = lancet_palette[9]),
        plot.subtitle = element_text(size = 7.2, hjust = 0, colour = lancet_palette[8]),
        axis.title.x = element_text(size = 8, colour = lancet_palette[9], margin = margin(t = 4)),
        axis.text.x = element_text(size = 7, colour = lancet_palette[9]),
        plot.margin = margin(6, 9, 6, 9))

out_dir <- file.path(paths$figures, "panel_A_corrected_workflow")
outputs <- export_figure5_plot(p, file.path(out_dir, "figure5_six_panel_A_corrected_workflow"), 13.5, 4.9)
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5_panel_A.rds"))
figure5_write_tsv(workflow[, .(step, title)], file.path(paths$metadata, "figure5_six_panel_A_workflow_nodes.tsv"))
figure5_write_tsv(data.table::data.table(step = seq_along(workflow_formal), workflow_node = workflow_formal),
                  file.path(paths$metadata, "figure5_six_panel_A_workflow_nodes_formal.tsv"))
figure5_write_tsv(state_labels, file.path(paths$metadata, "figure5_six_panel_A_state_anchors.tsv"))
figure5_write_json(list(
  panel = "5A",
  title = "Corrected temporal-positioning analysis workflow",
  workflow_nodes = workflow$title,
  workflow_nodes_formal = workflow_formal,
  state_anchors = state_labels$state,
  state_anchor_note = "State labels provide descriptive anchors and do not define discrete temporal transitions.",
  excluded_from_main_frame = c("independent trajectory-method rows", "independent resampling rows", "pre-drawn axis stage bands"),
  extended_data_audit_label = "Sensitivity and coverage audits—Extended Data",
  bootstrap_coverage_note = "Because eligible patient-token data did not span the complete pseudotemporal range, the primary bootstrap used the prespecified sample-token coverage fallback.",
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_A_report.json"))
