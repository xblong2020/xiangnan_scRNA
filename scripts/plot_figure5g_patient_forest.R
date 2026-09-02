#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
paths <- figure5_onset_fix_paths(root)
patients <- fread(file.path(paths$metadata, "figure5g_patient_temporal_differences.tsv"))
meta <- fread(file.path(paths$metadata, "figure5g_patient_meta_analysis.tsv"))[stratum == "overall" & status == "fitted"]
pooled <- if (nrow(meta) && "pooled_delta" %in% names(meta)) meta[, .(patient_id = "Random-effects pooled", dataset_id = "All", comparison, delta = pooled_delta, se = NA_real_, ci_lower, ci_upper,
                   n_cells = NA_integer_, n_bins = NA_integer_, same_direction = pooled_delta > 0, row_type = "pooled")] else data.table()
patients[, row_type := "patient"]
patients <- patients[is.finite(delta) & is.finite(ci_lower) & is.finite(ci_upper)]
plot_dt <- rbindlist(list(patients, pooled), fill = TRUE)
plot_dt[, comparison := factor(comparison, levels = c("A to B", "B to C", "A to C"))]
plot_dt[, label := paste0(patient_id, ifelse(row_type == "pooled", "", paste0("  (", dataset_id, ")")))]
plot_dt[, label := factor(label, levels = rev(unique(label)))]
comparison_palette <- c("A to B" = axis_palette[["identity_loss"]], "B to C" = axis_palette[["stress_transition"]], "A to C" = lancet_palette[5])
if (nrow(plot_dt)) {
  p <- ggplot(plot_dt, aes(delta, label, colour = comparison)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = lancet_palette[8]) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.18, linewidth = 0.45) +
    geom_point(aes(shape = row_type), size = 1.8) + facet_wrap(~comparison, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = comparison_palette, guide = "none") + scale_shape_manual(values = c(patient = 16, pooled = 18), guide = "none") +
    labs(title = "G  Evaluable patient-token-level temporal ordering", subtitle = "Delta > 0.005 indicates the upstream axis is earlier; intervals use leave-one-bin-out jackknife",
         x = "Pseudotime difference (Delta)", y = NULL) + theme_figure5(6.8) + theme(axis.text.y = element_text(size = 5.8))
} else {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.58, label = "Not available", fontface = "bold", size = 4, colour = lancet_palette[8]) +
    annotate("text", x = 0.5, y = 0.43,
             label = "No patient token met both prespecified criteria:\n≥5 pseudotime bins and coverage ≤0.10 to ≥0.90",
             size = 3, lineheight = 1.05, colour = lancet_palette[9]) +
    scale_x_continuous(limits = c(0, 1), breaks = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    labs(title = "G  Evaluable patient-token-level temporal ordering",
         subtitle = "Coverage-limited result; no independent clinical patient validation is claimed",
         x = NULL, y = NULL) +
    theme_figure5() + theme(axis.line = element_blank(), axis.ticks = element_blank())
}
out_dir <- file.path(paths$figures, "figure5g_patient_temporal_forest")
export_figure5_plot(p, file.path(out_dir, "figure5g_patient_temporal_forest"), 7.2, max(4.4, uniqueN(plot_dt$label) * 0.22))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5g.rds"))
