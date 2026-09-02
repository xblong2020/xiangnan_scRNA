#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_paths(root)
cells <- readRDS(file.path(paths$processed, "figure5_oriented_cell_scores.rds"))
cells <- as.data.table(cells)[is.finite(pseudotime)]
cells[, pseudotime_bin := pmin(10L, floor(pseudotime * 10) + 1L)]

score_cols <- grep("(_score$|_score_(tf_expression|regulon_auc|celloracle_target|intersection|no_cell_cycle|no_generic)$)", names(cells), value = TRUE)
score_cols <- unique(c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score",
                       "cellrank_fate_prob_cnv_supported_malignant", "cnv_score", "proliferation_score", score_cols))
score_cols <- score_cols[score_cols %in% names(cells)]

aggregate_bins <- function(dt, unit) {
  group <- if (unit == "patient") c("method", "dataset_id", "patient_id", "patient_id_source", "patient_meta_eligible", "pseudotime_bin") else
    c("method", "dataset_id", "patient_id", "patient_id_source", "patient_meta_eligible", "sample_id", "pseudotime_bin")
  out <- dt[, c(list(
    pseudotime = median(pseudotime, na.rm = TRUE), n_cells = uniqueN(cell_id), n_states_bin = uniqueN(cell_disease_stage),
    state_composition = paste(names(sort(table(cell_disease_stage), decreasing = TRUE)), collapse = ";"),
    cnv_strict_fraction = mean(cnv_strict %in% TRUE, na.rm = TRUE)
  ), lapply(.SD, function(x) median(as.numeric(x), na.rm = TRUE))), by = group, .SDcols = score_cols]
  out[, aggregation_unit := unit]
  out
}

patient <- aggregate_bins(cells, "patient")
sample <- aggregate_bins(cells, "sample")

eligibility <- cells[, .(n_cells_total = uniqueN(cell_id), n_bins_total = uniqueN(pseudotime_bin), n_states_total = uniqueN(cell_disease_stage),
                         patient_meta_eligible = all(patient_meta_eligible)), by = .(method, patient_id)]
eligibility[, eligible_patient := patient_eligibility(n_cells_total, n_bins_total, n_states_total, patient_meta_eligible)]
patient <- merge(patient, eligibility, by = c("method", "patient_id"), all.x = TRUE, suffixes = c("", ".elig"))
sample <- merge(sample, eligibility, by = c("method", "patient_id"), all.x = TRUE, suffixes = c("", ".elig"))
pseudo <- rbindlist(list(patient, sample), fill = TRUE, use.names = TRUE)
setcolorder(pseudo, c("aggregation_unit", "method", "dataset_id", "patient_id", "sample_id", "pseudotime_bin", "pseudotime", "n_cells", "eligible_patient", setdiff(names(pseudo), c("aggregation_unit", "method", "dataset_id", "patient_id", "sample_id", "pseudotime_bin", "pseudotime", "n_cells", "eligible_patient"))))
figure5_write_tsv(pseudo, file.path(paths$metadata, "figure5_patient_pseudotime_pseudobulk.tsv.gz"))
figure5_write_tsv(eligibility, file.path(paths$metadata, "figure5_patient_eligibility.tsv"))
saveRDS(pseudo, file.path(paths$processed, "figure5_patient_pseudotime_pseudobulk.rds"))
figure5_write_json(list(n_bins = 10, primary_unit = "patient", secondary_unit = "sample", eligibility = list(min_cells = 50, min_bins = 3, min_states = 2),
                        n_eligible_by_method = as.list(eligibility[eligible_patient == TRUE, uniqueN(patient_id), by = method][, setNames(V1, method)])),
                   file.path(paths$metadata, "figure5_patient_pseudobulk_report.json"))
message("Patient/sample pseudobulk complete: ", nrow(pseudo), " rows.")
