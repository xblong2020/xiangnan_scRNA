#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
features <- fread(file.path(source_paths$metadata, "figure5_heatmap_cell_feature_values.tsv.gz"), showProgress = FALSE)
oriented <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_oriented_cell_scores.rds")))[method == "main/consensus pseudotime",
  .(cell_id, pseudotime, dataset_id, patient_id, patient_meta_eligible, cell_disease_stage)]
dt <- merge(features, oriented, by = "cell_id")
dt <- dt[is.finite(value) & is.finite(pseudotime)]
dt[, pseudotime_bin := pmin(10L, floor(pseudotime * 10) + 1L)]
pseudo <- dt[, .(value = median(value, na.rm = TRUE), pseudotime = median(pseudotime), n_cells = uniqueN(cell_id),
                 n_states = uniqueN(cell_disease_stage)), by = .(entity, axis, entity_type, dataset_id, patient_id, patient_meta_eligible, pseudotime_bin)]
pseudo <- pseudo[patient_meta_eligible == TRUE]

entities <- unique(pseudo[, .(entity, axis, entity_type)])
pred_rows <- list(); landmark_rows <- list()
for (i in seq_len(nrow(entities))) {
  meta <- entities[i]
  sub <- pseudo[entity == meta$entity]
  fit <- landmarks_from_table(sub, "value", meta$axis, k = 5L, adjusted = FALSE, use_random_effect = TRUE)
  pred <- copy(fit$predictions)
  if (all(!is.finite(pred$fit))) next
  pred[, row_z := as.numeric(scale(fit))]
  pred[, `:=`(entity = meta$entity, axis = meta$axis, entity_type = meta$entity_type)]
  pred_rows[[meta$entity]] <- pred
  direction_fit <- pred$fit
  if (meta$axis == "identity_loss" && !grepl("Identity loss module", meta$entity, fixed = TRUE)) direction_fit <- -direction_fit
  raw_baseline <- sub[pseudotime >= figure5_temporal_parameters$baseline_start &
                        pseudotime <= figure5_temporal_parameters$baseline_end, value]
  baseline <- resolve_temporal_baseline_scale(raw_baseline, direction_fit)
  lm <- estimate_temporal_landmarks(
    pred$pseudotime, direction_fit, meta$axis,
    baseline_scale = baseline$scale, baseline_scale_source = baseline$source,
    observed_range = range(sub$pseudotime, na.rm = TRUE)
  )
  diag <- attr(lm, "diagnostics")
  landmark_rows[[meta$entity]] <- data.table(entity = meta$entity, axis = meta$axis, entity_type = meta$entity_type,
                                              onset_time = lm$onset_time, t10 = lm$t10, t50 = lm$t50,
                                              maximum_slope_time = lm$maximum_slope_time, extremum_time = lm$extremum_time,
                                              peak_time = lm$peak_time, plateau_time = lm$plateau_time, decline_onset = lm$decline_onset,
                                              baseline_scale = diag$baseline_scale,
                                              baseline_scale_source = diag$baseline_scale_source,
                                              coverage_ok = diag$coverage_ok, failure_reason = diag$failure_reason)
}
predictions <- rbindlist(pred_rows, fill = TRUE)
landmarks <- rbindlist(landmark_rows, fill = TRUE)
landmarks[, axis_order := match(axis, names(figure5_axis_score_columns))]
setorder(landmarks, axis_order, maximum_slope_time, entity)
landmarks[, display_order := seq_len(.N)]
predictions <- merge(predictions, landmarks[, .(entity, display_order)], by = "entity")

figure5_write_tsv(predictions[, .(entity, axis, entity_type, display_order, pseudotime, fitted_value = fit, row_z)],
                  file.path(paths$metadata, "figure5c_heatmap_matrix.tsv.gz"))
figure5_write_tsv(landmarks, file.path(paths$metadata, "figure5c_temporal_landmarks.tsv"))
manifest <- fread(file.path(source_paths$metadata, "figure5_heatmap_gene_manifest.tsv"))
figure5_write_tsv(manifest, file.path(paths$metadata, "figure5c_heatmap_gene_manifest.tsv"))
saveRDS(predictions, file.path(paths$processed, "figure5c_heatmap_predictions.rds"))
figure5_write_json(list(n_entities = nrow(landmarks), n_expression_genes = uniqueN(landmarks[entity_type == "TF expression" | entity_type == "target gene"]$entity),
                        model = "patient-pseudobulk GAM followed by row-wise z-score", ordering = "maximum-slope time; no pseudotime-based gene selection",
                        temporal_parameters = figure5_temporal_parameters,
                        frozen_manifest = file.path(source_paths$metadata, "figure5_frozen_signature_audit.tsv")),
                   file.path(paths$metadata, "figure5c_temporal_heatmap_report.json"))
message("Temporal heatmap data prepared for ", nrow(landmarks), " frozen entities.")
