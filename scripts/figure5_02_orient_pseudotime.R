#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_paths(root)
cells <- readRDS(file.path(paths$processed, "figure5_three_axis_cell_scores.rds"))
cells <- as.data.table(cells)

method_columns <- c(
  "main/consensus pseudotime" = "driver_main_strict__pseudotime_median",
  "Monocle3" = "driver_main_strict__monocle3_pseudotime",
  "Slingshot scanVI" = "driver_main_strict__slingshot_scanvi_pseudotime",
  "Slingshot hepatocyte PCA" = "driver_main_strict__slingshot_hepatocyte_pca_pseudotime",
  "CytoTRACE2" = "cytotrace2_score"
)

long_rows <- list()
audit_rows <- list()
for (method in names(method_columns)) {
  column <- method_columns[[method]]
  pt <- if (column %in% names(cells)) cells[[column]] else rep(NA_real_, nrow(cells))
  malignant <- as.numeric(as.logical(cells$driver_primary_module3_cnv_supported) | grepl("malignant", cells$trajectory_role, ignore.case = TRUE))
  orientation <- decide_pseudotime_orientation(pt, cells$identity_program_score_original, cells$cnv_score,
                                                cells$cellrank_fate_prob_cnv_supported_malignant, malignant)
  oriented <- orientation$pseudotime_oriented
  audit_rows[[method]] <- data.table(
    method = method, source_column = column, status = if (sum(is.finite(pt)) >= 50) "Available" else "Not available",
    n_cells = sum(is.finite(pt)), flipped = orientation$flipped, original_orientation_score = orientation$original_score,
    oriented_score = orientation$oriented_score,
    rho_identity_loss = unname(orientation$evidence[["identity_loss"]]), rho_cnv = unname(orientation$evidence[["cnv"]]),
    rho_malignant_fate = unname(orientation$evidence[["malignant_fate"]]), rho_malignant_fraction = unname(orientation$evidence[["malignant_fraction"]]),
    start_state = if (sum(is.finite(oriented))) cells[which.min(oriented)]$cell_disease_stage else "",
    end_state = if (sum(is.finite(oriented))) cells[which.max(oriented)]$cell_disease_stage else ""
  )
  if (sum(is.finite(oriented)) >= 50) {
    out <- copy(cells)
    out[, `:=`(method = method, pseudotime_original = figure5_scale01(pt), pseudotime = oriented)]
    long_rows[[method]] <- out
  }
}

audit <- rbindlist(audit_rows, fill = TRUE)
audit <- rbindlist(list(audit, data.table(method = c("DPT", "CellRank pseudotime"), source_column = "", status = "Not available",
                                          n_cells = 0L, flipped = NA, original_orientation_score = NA_real_, oriented_score = NA_real_,
                                          rho_identity_loss = NA_real_, rho_cnv = NA_real_, rho_malignant_fate = NA_real_,
                                          rho_malignant_fraction = NA_real_, start_state = "", end_state = "")), fill = TRUE)
oriented <- rbindlist(long_rows, fill = TRUE)
figure5_write_tsv(audit, file.path(paths$metadata, "figure5_pseudotime_orientation_audit.tsv"))
figure5_write_tsv(oriented, file.path(paths$metadata, "figure5_oriented_cell_scores.tsv.gz"))
saveRDS(oriented, file.path(paths$processed, "figure5_oriented_cell_scores.rds"))
figure5_write_json(list(definition = "0=normal/reference or mature end; 1=CNV-supported malignant/malignant-like end",
                        methods = split(audit, seq_len(nrow(audit))), unavailable_are_not_substituted = TRUE),
                   file.path(paths$metadata, "figure5_pseudotime_orientation_report.json"))
message("Pseudotime orientation audit complete for ", nrow(audit), " methods/status rows.")
