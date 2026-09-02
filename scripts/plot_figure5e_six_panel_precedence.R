#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_onset_fix_paths(root, create = FALSE)
old <- fread(file.path(source_paths$metadata, "figure5e_precedence_probabilities.tsv"))
old <- old[!is.na(landmark)]
boot <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))

# The corrected bootstrap did not previously export t10 precedence.  Derive it
# from the same resamples and tie rule used for onset, t50 and maximum slope.
make_rows <- function(landmark_label, bootstrap_landmark, onset_grade = FALSE) {
  wide <- dcast(boot[landmark == bootstrap_landmark], iteration ~ axis, value.var = "time")
  onset_fraction <- boot[landmark == "onset_time", .(fraction = mean(is.finite(time))), by = axis]
  onset_fraction <- setNames(onset_fraction$fraction, onset_fraction$axis)
  out <- lapply(seq_len(nrow(precedence_pairs)), function(i) {
    pair <- precedence_pairs[i]
    statistics <- tie_aware_precedence_probability(wide[[pair$upstream_axis]], wide[[pair$downstream_axis]])
    onset_pair_fraction <- min(onset_fraction[[pair$upstream_axis]], onset_fraction[[pair$downstream_axis]], na.rm = TRUE)
    grade <- classify_precedence_evidence(
      statistics$probability, statistics$valid_fraction, statistics$tie_fraction,
      statistics$delta_q025, statistics$delta_q975,
      onset_finite_fraction = if (onset_grade) onset_pair_fraction else NA_real_
    )
    mc_lower <- if (is.finite(statistics$Monte_Carlo_SE)) max(0, statistics$probability - 1.96 * statistics$Monte_Carlo_SE) else NA_real_
    mc_upper <- if (is.finite(statistics$Monte_Carlo_SE)) min(1, statistics$probability + 1.96 * statistics$Monte_Carlo_SE) else NA_real_
    data.table(
      landmark = landmark_label,
      comparison = pair$comparison,
      upstream_axis = pair$upstream_axis,
      downstream_axis = pair$downstream_axis,
      n_bootstrap_total = uniqueN(boot$iteration),
      n_valid = statistics$n_valid,
      n_earlier = statistics$n_earlier,
      n_tied = statistics$n_tied,
      n_later = statistics$n_later,
      valid_fraction = statistics$valid_fraction,
      tie_fraction = statistics$tie_fraction,
      probability = statistics$probability,
      median_delta = statistics$median_delta,
      delta_q025 = statistics$delta_q025,
      delta_q975 = statistics$delta_q975,
      Monte_Carlo_SE = statistics$Monte_Carlo_SE,
      mc_interval_lower = mc_lower,
      mc_interval_upper = mc_upper,
      mc_interval_label = "Monte Carlo interval only",
      onset_finite_fraction_pair = onset_pair_fraction,
      evidence_grade = grade,
      is_primary = landmark_label == "maximum_slope"
    )
  })
  rbindlist(out, fill = TRUE)
}

t10_rows <- make_rows("t10", "t10")
dt <- rbindlist(list(old, t10_rows), fill = TRUE)
dt <- dt[landmark %chin% c("onset", "t10", "t50", "maximum_slope")]
dt[, landmark_label := factor(landmark,
                              levels = c("onset", "t10", "t50", "maximum_slope"),
                              labels = c("Corrected onset", "t10", "t50", "Maximum slope"))]
dt[, comparison := factor(comparison, levels = c("A before B", "B before C", "A before C"))]
dt[, label := ifelse(is.finite(probability),
                     sprintf("Pr(earlier)=%.2f\nΔ[%.2f, %.2f]\ntie=%.1f%%",
                             probability, delta_q025, delta_q975, 100 * tie_fraction),
                     "Not available")]

p <- ggplot(dt, aes(comparison, landmark_label, fill = evidence_grade)) +
  geom_tile(colour = lancet_palette[8], linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.25, lineheight = 0.88) +
  scale_fill_manual(values = evidence_palette, drop = FALSE) +
  labs(title = "E  Tie-aware temporal precedence across complementary landmarks",
       subtitle = "Cells report Pr(earlier), bootstrap delta quantiles and tie fraction; peak is retained in D only",
       x = NULL, y = NULL, fill = "Evidence") +
  theme_figure5() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

out_dir <- file.path(paths$figures, "panel_E_precedence")
outputs <- export_figure5_plot(p, file.path(out_dir, "figure5_six_panel_E_precedence"), 6.2, 4.3)
figure5_write_tsv(dt, file.path(paths$metadata, "figure5_six_panel_E_precedence_probabilities.tsv"))
figure5_write_json(list(
  panel = "5E",
  title = "Tie-aware temporal precedence across complementary landmarks",
  probability_label = "Pr(earlier)",
  landmarks = c("corrected onset", "t10", "t50", "maximum slope"),
  peak_policy = "Peak is retained in Figure 5D only unless a prespecified interior stable peak audit passes; no such audit passed.",
  tie_tolerance = figure5_temporal_parameters$precedence_tolerance,
  source_namespace = "figure5_temporal_positioning_onset_fix",
  outputs = as.list(outputs)
), file.path(paths$metadata, "figure5_six_panel_E_report.json"))
dir.create(file.path(paths$processed, "plot_objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$processed, "plot_objects", "figure5_panel_E.rds"))

