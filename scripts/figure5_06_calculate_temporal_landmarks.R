#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
pseudo <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
selection <- select_figure5_primary_pseudobulk(pseudo)
primary <- selection$data
rows <- lapply(names(figure5_axis_score_columns), function(axis) {
  axis_name <- axis
  fit <- landmarks_from_table(primary, figure5_axis_score_columns[[axis_name]], axis_name,
                              k = 5L, adjusted = FALSE, use_random_effect = TRUE)
  lm <- fit$landmarks
  data.table(axis = axis_name, landmark = names(lm), time = as.numeric(unlist(lm)), method = "main/consensus pseudotime",
             model = "patient-pseudobulk GAM", scale = "0-1 oriented pseudotime")
})
landmarks <- rbindlist(rows)
diagnostics <- rbindlist(lapply(names(figure5_axis_score_columns), function(axis_name) {
  fit <- landmarks_from_table(primary, figure5_axis_score_columns[[axis_name]], axis_name,
                              k = 5L, adjusted = FALSE, use_random_effect = TRUE)
  as.data.table(fit$diagnostics)[, axis := axis_name]
}), fill = TRUE)
figure5_write_tsv(landmarks, file.path(paths$metadata, "figure5_temporal_landmarks.tsv"))
figure5_write_tsv(diagnostics, file.path(paths$metadata, "figure5_temporal_landmark_diagnostics.tsv"))
figure5_write_json(list(definitions = list(onset = "first post-baseline sustained positive run meeting effect and derivative thresholds",
                                            t10 = "first persistent 10% baseline-to-post-baseline-peak crossing",
                                            t50 = "first persistent 50% baseline-to-post-baseline-peak crossing",
                                            maximum_slope = "maximum positive central derivative after 0.10 and before min(peak,0.95)",
                                            extremum = "post-baseline fitted peak for all high-score programme axes",
                                            plateau = "first sustained low-absolute-derivative run after maximum slope",
                                            decline = "first sustained negative-derivative run after peak"),
                        fixed_parameters = figure5_temporal_parameters,
                        analysis_unit = selection$analysis_unit_note,
                        coverage_fallback_reason = selection$fallback_reason,
                        diagnostics = split(diagnostics, diagnostics$axis),
                        landmarks = split(landmarks, landmarks$axis)), file.path(paths$metadata, "figure5_temporal_landmarks_report.json"))
