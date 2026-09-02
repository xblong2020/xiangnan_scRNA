#!/usr/bin/env Rscript

## Figure 3A: overlapping stress-transition phases and empirical candidate matrix.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(scales)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "scripts", "figure3_egr1_common.R"))

data_dir <- normalizePath(
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3a_stress_transition")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3a_stress_transition")),
  mustWork = FALSE
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(data_dir, "figure3a_candidate_evidence_matrix.tsv")
if (!file.exists(matrix_path)) stop("Figure 3A candidate matrix is missing. Run prepare_figure3a_stress_transition_selection.py.")
candidate <- read.delim(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)

phases <- data.frame(
  xmin = c(0.05, 0.75, 1.55, 2.60),
  xmax = c(1.55, 2.40, 3.25, 4.05),
  ymin = c(3.25, 2.30, 1.35, 0.40),
  ymax = c(3.90, 2.95, 2.00, 1.05),
  fill = c("#B8B8B8", "#56B4E9", "#E69F00", "#D55E00"),
  label = c(
    "HNF4A/PPARA-associated\nhepatocyte identity attenuation",
    "AP-1 rapid stress response\nJUN · JUNB · JUND · FOS · FOSB · ATF3",
    "CEBPB/EGR1-associated\ntransition programme",
    "SOX4-associated\nmalignant-state stabilization"
  ),
  stringsAsFactors = FALSE
)

connectors <- data.frame(
  x = c(1.15, 1.90, 2.85),
  y = c(3.23, 2.28, 1.33),
  xend = c(1.45, 2.20, 3.10),
  yend = c(2.98, 2.03, 1.08)
)

p_architecture <- ggplot() +
  geom_rect(
    data = phases,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.58,
    colour = "white",
    linewidth = 0.7
  ) +
  geom_text(
    data = phases,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    family = "sans",
    size = 2.45,
    lineheight = 0.92,
    colour = "grey12"
  ) +
  geom_curve(
    data = connectors,
    aes(x = x, y = y, xend = xend, yend = yend),
    curvature = 0.25,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.45,
    arrow = arrow(length = grid::unit(0.045, "inches"), type = "open")
  ) +
  annotate(
    "label",
    x = 3.42, y = 1.50,
    label = "EGR1\nprincipal perturbation representative",
    size = 2.15,
    family = "sans",
    linewidth = 0.35,
    fill = "#FFF5CC",
    colour = "#7A4E00"
  ) +
  annotate(
    "text",
    x = 2.02, y = 4.27,
    label = "overlapping regulatory phases",
    family = "sans",
    fontface = "italic",
    size = 2.7,
    colour = "grey25"
  ) +
  annotate(
    "text",
    x = 2.02, y = 0.08,
    label = "Partially ordered architecture; dashed connectors indicate phase overlap, not a proven linear cascade.",
    family = "sans",
    size = 2.25,
    colour = "grey30"
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.05, 4.15), ylim = c(-0.05, 4.45), clip = "off") +
  labs(title = "AP-1/CEBPB/EGR1 stress-transition programme") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 6)),
    plot.margin = margin(6, 4, 6, 4)
  )

metric_columns <- c(
  celloracle_evidence = "CellOracle",
  sctenifoldknk_evidence = "scTenifoldKnk",
  cross_method_concordance = "Concordance",
  transition_state_specificity = "Transition\nspecificity",
  temporal_positioning = "Temporal\nposition",
  leave_one_dataset_out_stability = "LODO/LOSO\nstability",
  proliferation_dependency = "Proliferation\ndep.",
  generic_stress_risk = "Generic stress\nrisk",
  literature_overlap = "Literature\noverlap",
  selection_score_scaled = "Selection\nscore"
)

long_rows <- lapply(names(metric_columns), function(column) {
  data.frame(
    candidate = candidate$candidate,
    metric = unname(metric_columns[[column]]),
    score = as.numeric(candidate[[column]]),
    label = sprintf("%.2f", as.numeric(candidate[[column]])),
    stringsAsFactors = FALSE
  )
})
long <- do.call(rbind, long_rows)
candidate_order <- candidate$candidate[order(candidate$selection_rank)]
long$candidate <- factor(long$candidate, levels = rev(candidate_order))
long$metric <- factor(long$metric, levels = unname(metric_columns))
roles <- candidate[, c("candidate", "final_role"), drop = FALSE]
roles$candidate <- factor(roles$candidate, levels = rev(candidate_order))

p_heatmap <- ggplot(long, aes(x = metric, y = candidate, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = label), size = 2.25, family = "sans", colour = "grey10") +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#56B4E9", "#009E73"),
    limits = c(0, 1),
    oob = scales::squish,
    name = "Standardized\nmetric"
  ) +
  labs(x = NULL, y = NULL, title = "Candidate evidence matrix") +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 6)),
    axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 7.5, colour = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.42, "cm"),
    plot.margin = margin(6, 2, 6, 4)
  )

p_roles <- ggplot(roles, aes(x = 0, y = candidate, label = final_role)) +
  geom_text(hjust = 0, size = 2.35, family = "sans", lineheight = 0.92) +
  xlim(0, 1) +
  labs(title = "Final role") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 8.5, hjust = 0, margin = margin(b = 6)),
    plot.margin = margin(6, 2, 6, 0)
  )

p_matrix <- p_heatmap + p_roles + plot_layout(widths = c(4.3, 1.25))
selection_sentence <- "EGR1 was selected as the principal perturbation representative of the AP-1/CEBPB/EGR1 stress-transition programme."
combined <- (p_architecture | p_matrix) +
  plot_layout(widths = c(0.88, 1.62)) +
  plot_annotation(
    caption = selection_sentence,
    theme = theme(
      plot.caption = element_text(size = 8, hjust = 0.5, family = "sans", colour = "grey15", margin = margin(t = 6))
    )
  )

stem <- "figure3a_stress_transition_selection"
figure3_save(combined, figure_dir, stem, width = 13.2, height = 5.3, tiff = TRUE)

report <- list(
  module = "Figure 3A plotting",
  target_tf = "EGR1",
  plotting_language = "R",
  r_version = R.version.string,
  input = figure3_norm_path(matrix_path),
  n_candidates = nrow(candidate),
  empirical_ranking = as.list(candidate$candidate[order(candidate$selection_rank)]),
  architecture = "overlapping regulatory phases with dashed connectors; not a strict causal cascade",
  ap1_members_shown = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "ATF3"),
  outputs = list(
    pdf = figure3_norm_path(file.path(figure_dir, paste0(stem, ".pdf"))),
    png = figure3_norm_path(file.path(figure_dir, paste0(stem, ".png"))),
    svg = figure3_norm_path(file.path(figure_dir, paste0(stem, ".svg"))),
    tiff = figure3_norm_path(file.path(figure_dir, paste0(stem, ".tiff")))
  ),
  caveat = "The evidence matrix documents perturbation-target selection. It does not establish EGR1 as a proven causal driver, a unique effective gene, or part of an obligatory linear cascade."
)
figure3_write_json(report, file.path(data_dir, "figure3a_r_plot_report.json"))
message("Figure 3A written to: ", figure_dir)
