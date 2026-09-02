#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

effects <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_effects.tsv.gz"))
display_rows <- c("HNF4A KO", "HNF4A restore/OE", "PPARA KO", "PPARA restore/OE", "EGR1 KO", "CEBPB KO",
  "AP-1 member-KO aggregate", "SOX4 KO", "SOX4 OE", "HLF KO", "IRF1 KO", "MAFB KO", "MAFF KO", "MYC KO")
plot_data <- effects[perturbation %in% display_rows]
plot_data[, perturbation := factor(perturbation, levels = rev(display_rows))]
plot_data[, output := factor(output, levels = FIGURE6_CORE_OUTPUTS, labels = unname(FIGURE6_OUTPUT_LABELS[FIGURE6_CORE_OUTPUTS]))]
plot_data[, sig := figure6_significance(fdr)]
plot_data[, unstable := availability == "Available" & is.finite(stability) & stability < 0.60]
plot_data[, alpha_value := fifelse(availability != "Available", 1, pmax(0.45, stability))]
plot_data[, axis_label := factor(axis, levels = c("identity_axis", "stress_axis", "sox4_axis", "control"),
  labels = c("Axis A", "Axis B", "Axis C", "Control"))]
figure6_fwrite(plot_data, file.path(FIGURE6_METADATA_DIR, "figure6a_matrix_plot_data.tsv"), compress = FALSE)

limit <- max(abs(plot_data$effect_estimate), na.rm = TRUE)
limit <- ifelse(is.finite(limit) && limit > 0, limit, 1)
p <- ggplot(plot_data, aes(output, perturbation)) +
  geom_tile(aes(fill = effect_estimate, alpha = alpha_value), colour = "white", linewidth = 0.55) +
  geom_text(aes(label = fifelse(availability != "Available", "×", sig)),
    colour = dark_text, size = 3.0, fontface = "bold") +
  geom_point(data = plot_data[unstable == TRUE], shape = 4, size = 1.4, stroke = 0.35, colour = neutral_gray) +
  scale_fill_gradient2(low = effect_gradient["low"], mid = effect_gradient["mid"], high = effect_gradient["high"],
    midpoint = 0, limits = c(-limit, limit), oob = scales::squish, na.value = "white", name = "Standardized\neffect") +
  scale_alpha_identity() +
  facet_grid(axis_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Perturbation-response matrix", x = NULL, y = NULL,
    caption = "Stars: sample-level dataset-stratified bootstrap FDR. ×: Not available. Small × overlay: dataset stability <0.60. AP-1 is a member-KO median aggregate.") +
  figure6_theme() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1), strip.placement = "outside",
    strip.background = element_blank(), strip.text.y.left = element_text(angle = 90, size = 7, colour = dark_text),
    panel.spacing.y = grid::unit(0.08, "cm"), legend.position = "right")

out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6a_perturbation_matrix")
figure6_save(p, out_dir, "figure6a_perturbation_response_matrix", 8.3, 6.7)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6a_plot.rds"))
figure6_write_json(list(
  panel = "Figure 6A", plot = "ggplot2::geom_tile", midpoint = 0, n_rows = nrow(plot_data),
  unavailable_cells = sum(plot_data$availability != "Available"), fdr_source = "sample-level bootstrap"
), file.path(FIGURE6_METADATA_DIR, "figure6a_matrix_report.json"))

