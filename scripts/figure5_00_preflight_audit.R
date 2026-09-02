#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_paths(root)

inputs <- c(
  module9_cell_scores = file.path(root, "metadata/driver/module9_1_temporal_cell_scores.tsv.gz"),
  driver_cells = file.path(root, "metadata/driver/driver_module6_1_cells.tsv.gz"),
  trajectory_main = file.path(root, "metadata/trajectory/trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz"),
  cellrank_fate = file.path(root, "metadata/driver/driver_module6_2_cellrank_fate_probabilities.tsv.gz"),
  cistarget_auc = file.path(root, "metadata/driver/driver_module6_3c_cistarget_regulon_auc.tsv.gz"),
  pyscenic_auc = file.path(root, "metadata/driver/driver_module6_3_pyscenic_regulon_auc.tsv.gz"),
  celloracle_links = file.path(root, "metadata/driver/celloracle_module6_7_grn_links_filtered.tsv.gz"),
  sctenifold_targets = file.path(root, "metadata/driver/module8_tf_target_signature_genes.tsv"),
  cytotrace2 = file.path(root, "metadata/figure1c/figure1c_cytotrace2_scores_by_cell.hepatocyte.tsv.gz"),
  driver_h5ad = file.path(root, "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad")
)

file_rows <- rbindlist(lapply(names(inputs), function(name) {
  info <- file.info(inputs[[name]])
  data.table(category = "input_file", item = name, status = if (file.exists(inputs[[name]])) "Available" else "Not available",
             value = if (file.exists(inputs[[name]])) as.character(info$size) else "", details = inputs[[name]])
}))

required <- inputs[c("module9_cell_scores", "driver_cells", "trajectory_main", "cellrank_fate", "cistarget_auc", "sctenifold_targets", "driver_h5ad")]
if (any(!file.exists(required))) stop("Required frozen inputs are missing: ", paste(names(required)[!file.exists(required)], collapse = ", "))

message("Reading frozen Figure 1-4/module9 inputs for exact cell audit...")
scores <- fread(inputs[["module9_cell_scores"]], showProgress = FALSE)
driver <- fread(inputs[["driver_cells"]], showProgress = FALSE)
fate <- fread(inputs[["cellrank_fate"]], showProgress = FALSE)
cyto <- if (file.exists(inputs[["cytotrace2"]])) fread(inputs[["cytotrace2"]], select = c("cell_id", "CytoTRACE2_Score"), showProgress = FALSE) else data.table()
strict <- unique(scores[run_id == "main_strict", .(cell_id, dataset, sample_id, study_sample, cnv_sample,
                                                        cell_disease_stage, trajectory_role)])
mapped <- derive_figure5_patient_id(strict$dataset, strict$cnv_sample)
strict <- cbind(strict, mapped)

counts <- rbindlist(list(
  data.table(category = "count", item = c("strict_cells", "datasets", "study_samples", "cnv_samples", "patient_tokens", "meta_eligible_patient_tokens"),
             status = "Available", value = as.character(c(nrow(strict), uniqueN(strict$dataset), uniqueN(strict$study_sample),
                                                           uniqueN(strict$cnv_sample), uniqueN(strict$patient_id),
                                                           uniqueN(strict[patient_meta_eligible == TRUE]$patient_id))), details = "main_strict frozen set"),
  strict[, .(category = "dataset_count", item = dataset, status = "Available", value = as.character(.N), details = "strict cells"), by = dataset][, dataset := NULL],
  strict[, .(category = "state_count", item = cell_disease_stage, status = "Available", value = as.character(.N), details = "strict cells"), by = cell_disease_stage][, cell_disease_stage := NULL]
), fill = TRUE)

match_rows <- data.table(
  category = "cell_id_match",
  item = c("driver_cells", "cellrank_fate", "cytotrace2"),
  status = c(if (all(strict$cell_id %in% driver$cell_id)) "Complete" else "Partial",
             if (all(strict$cell_id %in% fate$cell_id)) "Complete" else "Partial",
             if (nrow(cyto) && all(strict$cell_id %in% cyto$cell_id)) "Complete" else if (nrow(cyto)) "Partial" else "Not available"),
  value = as.character(c(sum(strict$cell_id %in% driver$cell_id), sum(strict$cell_id %in% fate$cell_id), sum(strict$cell_id %in% cyto$cell_id))),
  details = paste0("of ", nrow(strict), " strict cells")
)

methods <- data.table(
  category = "trajectory_method",
  item = c("main/consensus pseudotime", "DPT", "Monocle3", "Slingshot scanVI", "Slingshot hepatocyte PCA", "CellRank pseudotime", "CytoTRACE2", "CellRank malignant fate"),
  status = c("Available", "Not available", "Available", "Available", "Available", "Not available",
             if (nrow(cyto)) "Available" else "Not available", "Available"),
  value = "", details = c("median of three frozen methods", "no frozen DPT result found", "module5_3", "module5_3", "module5_3",
                           "fate probability exists; no CellRank pseudotime", "Figure 1C frozen score", "module6_2")
)

covariates <- data.table(
  category = "covariate", item = c("proliferation", "S_phase", "G2M", "hypoxia", "mitochondrial", "ribosomal", "dissociation_stress", "CNV_score", "strict_CNV_label"),
  status = c(if ("proliferation_score_z" %in% names(driver)) "Available" else "Not available",
             rep("Not available", 6), if ("cnv_proxy_z" %in% names(driver)) "Available" else "Not available",
             if ("driver_primary_module3_cnv_supported" %in% names(driver)) "Available" else "Not available"),
  value = "", details = "absence is retained as Not available"
)

patient_counts <- strict[, .(n_cells = .N, n_states = uniqueN(cell_disease_stage)), by = .(dataset, patient_id, patient_id_source, patient_meta_eligible)]
patient_rows <- patient_counts[, .(category = "patient_feasibility", item = patient_id,
                                   status = ifelse(n_cells >= 50 & n_states >= 2 & patient_meta_eligible, "Potentially eligible", "Insufficient/aggregate"),
                                   value = as.character(n_cells), details = paste0("dataset=", dataset, "; states=", n_states, "; source=", patient_id_source))]

preflight <- rbindlist(list(file_rows, counts, match_rows, methods, covariates, patient_rows), fill = TRUE, use.names = TRUE)
figure5_write_tsv(preflight, file.path(paths$metadata, "figure5_preflight_report.tsv"))

# Git cannot protect this untracked analysis tree, so record a byte-level Figure 1-4 baseline.
candidate_roots <- c(file.path(root, "scripts"), file.path(root, "metadata"), file.path(root, "figures"), file.path(root, "reports"))
all_files <- unlist(lapply(candidate_roots[file.exists(candidate_roots)], function(x) list.files(x, recursive = TRUE, full.names = TRUE)), use.names = FALSE)
rel <- substring(normalizePath(all_files, winslash = "/", mustWork = FALSE), nchar(normalizePath(root, winslash = "/")) + 2L)
keep <- grepl("figure[1-4]|fig[1-4]", rel, ignore.case = TRUE) & !grepl("figure5", rel, ignore.case = TRUE)
baseline_files <- all_files[keep & file.exists(all_files)]
baseline <- data.table(file_path = normalizePath(baseline_files, winslash = "/"), relative_path = rel[keep & file.exists(all_files)])
if (nrow(baseline)) {
  info <- file.info(baseline$file_path)
  baseline[, `:=`(size_bytes = info$size, modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE), md5 = unname(tools::md5sum(file_path)))]
}
figure5_write_tsv(baseline, file.path(paths$metadata, "figure5_figure1_4_baseline.tsv"))
saveRDS(list(scores = scores, driver = driver, fate = fate, cytotrace = cyto, strict = strict), file.path(paths$processed, "figure5_preflight_cache.rds"))

report <- list(
  module = "Figure 5 preflight", created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = as.list(inputs), n_strict_cells = nrow(strict), n_datasets = uniqueN(strict$dataset),
  n_study_samples = uniqueN(strict$study_sample), n_cnv_samples = uniqueN(strict$cnv_sample),
  n_patient_tokens = uniqueN(strict$patient_id), patient_id_policy = "explicit sample-name token; aggregate objects excluded from patient meta-analysis",
  available_methods = methods[status == "Available"]$item, unavailable_methods = methods[status == "Not available"]$item,
  figure1_4_baseline_files = nrow(baseline), review_risk_flags = c("patient_id is derived because no frozen patient_id field exists", "DPT and CellRank pseudotime are unavailable")
)
figure5_write_json(report, file.path(paths$metadata, "figure5_preflight_report.json"))
message("Figure 5 preflight complete: ", nrow(strict), " strict cells.")
