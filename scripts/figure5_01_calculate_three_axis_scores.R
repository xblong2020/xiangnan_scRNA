#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_paths(root)
set.seed(20260805)

cache_path <- file.path(paths$processed, "figure5_preflight_cache.rds")
if (!file.exists(cache_path)) stop("Run figure5_00_preflight_audit.R first")
cache <- readRDS(cache_path)
scores_long <- as.data.table(cache$scores)
base <- unique(scores_long[run_id == "main_strict" & method == "monocle3", .(
  cell_id, identity_program_score_original = A_hnf4a_ppara_retention,
  identity_loss_score_raw = A_hnf4a_ppara_loss,
  stress_transition_score_raw = B_transition_activation,
  sox4_stabilization_score_raw = C_sox4_axis,
  stress_ap1_score_raw = B_ap1_activation,
  stress_cebpb_egr1_score_raw = B_cebpb_egr1_activation,
  malignant_like_fate_score_raw = C_malignant_like_fate,
  cellrank_fate_prob_cnv_supported_malignant
)], by = "cell_id")

driver <- as.data.table(cache$driver)
driver_keep <- intersect(c("cell_id", "dataset", "sample_id", "study_sample", "cnv_sample", "sample_source_class",
                           "cell_disease_stage", "trajectory_role", "trajectory_root_end_role", "malignant_hcc_call",
                           "cnv_proxy_z", "cnv_proxy_burden", "hcc_malignant_associated_score_z", "proliferation_score_z",
                           "driver_main_strict__module_Proliferation",
                           "driver_primary_module3_cnv_supported", "driver_primary_cnv_evidence_tier",
                           "driver_main_strict__monocle3_pseudotime", "driver_main_strict__slingshot_scanvi_pseudotime",
                           "driver_main_strict__slingshot_hepatocyte_pca_pseudotime", "driver_main_strict__pseudotime_median",
                           "driver_main_strict__eligible", "trajectory_include_cnv_strict"), names(driver))
base <- merge(base, unique(driver[, ..driver_keep], by = "cell_id"), by = "cell_id", all.x = TRUE)
setnames(base, "dataset", "dataset_id")
base[, input_sample_id := sample_id]
base[, sample_id := paste(dataset_id, cnv_sample, sep = "::")]
base[, study_sample_id := study_sample]
mapped <- derive_figure5_patient_id(base$dataset_id, base$cnv_sample)
base[, `:=`(patient_id = mapped$patient_id, patient_id_source = mapped$patient_id_source,
            patient_meta_eligible = mapped$patient_meta_eligible)]
base[, `:=`(
  identity_loss_score = robust_z_by_group(identity_loss_score_raw, dataset_id),
  stress_transition_score = robust_z_by_group(stress_transition_score_raw, dataset_id),
  sox4_stabilization_score = robust_z_by_group(sox4_stabilization_score_raw, dataset_id),
  proliferation_score = fifelse(is.finite(proliferation_score_z), proliferation_score_z,
                                 robust_z_by_group(driver_main_strict__module_Proliferation, dataset_id)),
  cnv_score = fifelse(is.finite(cnv_proxy_z), cnv_proxy_z, hcc_malignant_associated_score_z),
  cnv_strict = as.logical(driver_primary_module3_cnv_supported) | trajectory_role == "normal_reference"
)]

if (nrow(cache$cytotrace)) {
  cyto <- unique(as.data.table(cache$cytotrace), by = "cell_id")
  setnames(cyto, "CytoTRACE2_Score", "cytotrace2_score")
  base <- merge(base, cyto, by = "cell_id", all.x = TRUE)
} else base[, cytotrace2_score := NA_real_]

# Read the frozen H5AD from R through reticulate; Python performs object I/O only and no plotting.
suppressPackageStartupMessages({library(reticulate); library(Matrix); library(UCell)})
python <- "C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe"
reticulate::use_python(python, required = TRUE)
ad <- reticulate::import("anndata", convert = FALSE)
h5ad_path <- normalizePath(file.path(root, "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad"), winslash = "/")
adata <- ad$read_h5ad(h5ad_path, backed = "r")
memory <- adata$to_memory()
x <- reticulate::py_to_r(memory$X)
obs_names <- reticulate::py_to_r(memory$obs_names$to_list())
var_names <- toupper(reticulate::py_to_r(memory$var_names$to_list()))
adata$file$close()
idx <- match(base$cell_id, obs_names)
if (anyNA(idx)) stop("Strict cells are missing from frozen driver H5AD: ", sum(is.na(idx)))
mat <- Matrix::t(x[idx, , drop = FALSE])
mat <- methods::as(mat, "CsparseMatrix")
rownames(mat) <- var_names
colnames(mat) <- base$cell_id

gene_sets <- readRDS(file.path(paths$processed, "figure5_frozen_gene_sets.rds"))
ucell_cache <- file.path(paths$processed, "figure5_ucell_sensitivity_scores.rds")
if (file.exists(ucell_cache)) {
  ucell <- readRDS(ucell_cache)
} else {
  ucell <- UCell::ScoreSignatures_UCell(mat, features = gene_sets, maxRank = min(1500L, nrow(mat)),
                                       missing_genes = "skip", ncores = 1L)
  saveRDS(ucell, ucell_cache)
}
ucell <- as.data.table(ucell, keep.rownames = "cell_id")
if (anyDuplicated(ucell$cell_id)) stop("UCell output contains duplicated cell IDs")
base <- merge(base, ucell, by = "cell_id", all.x = TRUE)

u <- function(name) paste0(name, "_UCell")
base[, `:=`(
  identity_loss_score_tf_expression = -robust_z_by_group(get(u("tf_identity_retention")), dataset_id),
  stress_transition_score_tf_expression = robust_z_by_group(get(u("tf_stress_transition")), dataset_id),
  sox4_stabilization_score_tf_expression = robust_z_by_group(get(u("tf_sox4_stabilization")), dataset_id),
  identity_loss_score_celloracle_target = -robust_z_by_group(get(u("celloracle_identity_retention")), dataset_id),
  stress_transition_score_celloracle_target = robust_z_by_group(get(u("celloracle_stress_transition")), dataset_id),
  sox4_stabilization_score_celloracle_target = robust_z_by_group(get(u("celloracle_sox4_stabilization")), dataset_id),
  identity_loss_score_intersection = -robust_z_by_group(get(u("intersection_identity_retention")), dataset_id),
  stress_transition_score_intersection = robust_z_by_group(get(u("intersection_stress_transition")), dataset_id),
  sox4_stabilization_score_intersection = robust_z_by_group(get(u("intersection_sox4_stabilization")), dataset_id),
  identity_loss_score_no_cell_cycle = -robust_z_by_group(get(u("identity_no_cell_cycle")), dataset_id),
  stress_transition_score_no_cell_cycle = robust_z_by_group(get(u("stress_no_cell_cycle")), dataset_id),
  sox4_stabilization_score_no_cell_cycle = robust_z_by_group(get(u("sox4_no_cell_cycle")), dataset_id),
  stress_transition_score_no_generic = robust_z_by_group(get(u("stress_no_generic")), dataset_id)
)]

# Regulon-only sensitivity scores use the frozen cisTarget AUC table.
auc <- fread(file.path(root, "metadata/driver/driver_module6_3c_cistarget_regulon_auc.tsv.gz"), showProgress = FALSE)
setkey(auc, cell_id)
auc <- auc[base$cell_id]
reg_mean <- function(tfs) {
  cols <- intersect(paste0(tfs, "(+)"), names(auc))
  if (!length(cols)) return(rep(NA_real_, nrow(base)))
  rowMeans(as.matrix(auc[, ..cols]), na.rm = TRUE)
}
base[, `:=`(
  identity_loss_score_regulon_auc = -robust_z_by_group(reg_mean(c("HNF4A", "PPARA")), dataset_id),
  stress_transition_score_regulon_auc = robust_z_by_group(reg_mean(c("JUN", "JUNB", "JUND", "FOS", "ATF3", "CEBPB", "EGR1")), dataset_id),
  sox4_stabilization_score_regulon_auc = robust_z_by_group(reg_mean("SOX4"), dataset_id)
)]

score_path <- file.path(paths$metadata, "figure5_three_axis_cell_scores.tsv.gz")
figure5_write_tsv(base, score_path)
saveRDS(base, file.path(paths$processed, "figure5_three_axis_cell_scores.rds"))

# Feature values for the frozen heatmap: expression, regulon AUC, and module scores.
manifest <- fread(file.path(paths$metadata, "figure5_heatmap_gene_manifest_preavailability.tsv"))
available_genes <- intersect(unique(manifest$gene), rownames(mat))
manifest[, available := gene %chin% available_genes]
figure5_write_tsv(manifest, file.path(paths$metadata, "figure5_heatmap_gene_manifest.tsv"))
expr <- as.data.table(as.matrix(Matrix::t(mat[available_genes, , drop = FALSE])), keep.rownames = "cell_id")
expr_long <- melt(expr, id.vars = "cell_id", variable.name = "entity", value.name = "value")
manifest_expr <- unique(manifest[available == TRUE, .(entity = gene, axis, entity_type)])
manifest_expr[, type_priority := fcase(entity_type == "TF expression", 1L, entity_type == "target gene", 2L, default = 3L)]
setorder(manifest_expr, entity, type_priority)
manifest_expr <- unique(manifest_expr, by = "entity")
expr_long <- merge(expr_long, manifest_expr[, .(entity, axis, entity_type)], by = "entity")

reg_tfs <- unique(unlist(list(c("HNF4A", "PPARA"), c("JUN", "JUNB", "JUND", "FOS", "ATF3", "CEBPB", "EGR1"), "SOX4")))
reg_cols <- intersect(paste0(reg_tfs, "(+)"), names(auc))
reg <- cbind(data.table(cell_id = base$cell_id), auc[, ..reg_cols])
reg_long <- melt(reg, id.vars = "cell_id", variable.name = "entity", value.name = "value")
reg_long[, tf := sub("\\(\\+\\)$", "", entity)]
reg_long[, axis := fcase(tf %chin% c("HNF4A", "PPARA"), "identity_loss", tf == "SOX4", "sox4_stabilization", default = "stress_transition")]
reg_long[, `:=`(entity = paste0(tf, " regulon"), entity_type = "regulon", tf = NULL)]
modules <- melt(base[, .(cell_id, identity_loss_score, stress_transition_score, sox4_stabilization_score)], id.vars = "cell_id", variable.name = "axis", value.name = "value")
modules[, axis := sub("_score$", "", axis)]
modules[, `:=`(entity = paste0(figure5_axis_labels[axis], " module"), entity_type = "module score")]
feature_long <- rbindlist(list(expr_long, reg_long[, .(cell_id, entity, value, axis, entity_type)], modules[, .(cell_id, entity, value, axis, entity_type)]), fill = TRUE)
figure5_write_tsv(feature_long, file.path(paths$metadata, "figure5_heatmap_cell_feature_values.tsv.gz"))

report <- list(method = "Frozen primary scores plus UCell sensitivity scoring", ucell_version = as.character(packageVersion("UCell")),
               h5ad = h5ad_path, n_cells = nrow(base), n_genes_ranked = nrow(mat), n_heatmap_genes_available = length(available_genes),
               identity_loss_definition = "-z(identity_program_score_original), preserving the original identity score",
               standardization = "dataset-wise median/MAD robust z-score, clipped to [-5,5]", random_seed = 20260805)
figure5_write_json(report, file.path(paths$metadata, "figure5_three_axis_scores_report.json"))
message("Three-axis score table written: ", nrow(base), " cells.")
