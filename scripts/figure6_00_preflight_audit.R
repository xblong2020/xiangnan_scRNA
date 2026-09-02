#!/usr/bin/env Rscript

## Figure 6 frozen-input and perturbation availability audit.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(), mustWork = FALSE)
source(file.path(PROJECT_ROOT, "scripts", "figure6_common.R"))

metadata_root <- file.path(PROJECT_ROOT, "metadata", "driver")
trajectory_root <- file.path(PROJECT_ROOT, "metadata", "trajectory")

paths <- list(
  celloracle_report = file.path(metadata_root, "celloracle_module6_8_perturbation_report.json"),
  celloracle_cells = file.path(metadata_root, "celloracle_module6_8_cell_shift_summary.tsv.gz"),
  celloracle_gene_delta = file.path(metadata_root, "celloracle_module6_8_top_gene_delta_by_state.tsv.gz"),
  celloracle_tf_delta = file.path(metadata_root, "celloracle_module6_8_tf_delta_summary.tsv"),
  celloracle_grid = file.path(metadata_root, "celloracle_module6_8_grid_arrows.tsv.gz"),
  celloracle_grn = file.path(metadata_root, "celloracle_module6_7_grn_links_filtered.tsv.gz"),
  celloracle_metadata = file.path(metadata_root, "celloracle_module6_6_cell_metadata.tsv.gz"),
  oracle_object = file.path(PROJECT_ROOT, "data", "processed", "driver", "celloracle_module6_7", "celloracle_module6_7_fitted.celloracle.oracle"),
  sct_main = file.path(metadata_root, "sctenifoldknk_module7_2_main_strict_perturbation_genes.tsv"),
  sct_contract = file.path(metadata_root, "sctenifoldknk_module7_2_main_strict_result_contract.tsv"),
  scenic_regulons = file.path(metadata_root, "driver_module6_3_pyscenic_regulons.tsv"),
  cistarget_regulons = file.path(metadata_root, "driver_module6_3c_cistarget_regulon_summary.tsv"),
  cellrank_fate = file.path(metadata_root, "driver_module6_2_cellrank_fate_probabilities.tsv.gz"),
  strict_cnv = file.path(metadata_root, "driver_module6_1_cells.tsv.gz"),
  figure5_temporal = file.path(metadata_root, "module9_1_onset_times.tsv"),
  figure5_order = file.path(metadata_root, "module9_1_bootstrap_order_tests.tsv"),
  module9_direction = file.path(metadata_root, "module9_2_network_directionality_matrix.tsv"),
  module9_restore = file.path(metadata_root, "module9_2_restore_availability.tsv"),
  module9_pseudobulk = file.path(metadata_root, "module9_3_scrna_pseudobulk_axis_scores.tsv.gz"),
  external_validation = file.path(metadata_root, "module8_external_validation_report.json"),
  frozen_module_genes = file.path(trajectory_root, "trajectory_module5_4_module_gene_availability.tsv"),
  hnf4a_vector = file.path(metadata_root, "figure2c_hnf4a", "figure2c_hnf4a_grid_umap.tsv.gz"),
  egr1_vector = file.path(metadata_root, "figure3c_egr1", "figure3c_egr1_grid_umap.tsv.gz"),
  sox4_vector = file.path(metadata_root, "figure2c_sox4", "figure2c_sox4_grid_umap.tsv.gz"),
  three_axis_audit = file.path(metadata_root, "three_axis_figure_consistency", "figure2_figure3_figure4_consistency_report.json")
)

report6 <- if (file.exists(paths$celloracle_report)) read_json(paths$celloracle_report, simplifyVector = TRUE) else list()
available_ko <- as.character(report6$result$perturbed_tfs %||% character())
if (!length(available_ko) && file.exists(paths$celloracle_tf_delta)) {
  available_ko <- unique(figure6_fread(paths$celloracle_tf_delta, select = "tf")$tf)
}
sct_contract <- figure6_fread(paths$sct_contract)
sct_tfs <- if ("tf" %in% names(sct_contract)) as.character(sct_contract$tf) else character()
cell_meta <- figure6_fread(paths$celloracle_metadata)
n_cells <- if (nrow(cell_meta)) nrow(cell_meta) else as.integer(report6$result$n_cells %||% NA)
n_datasets <- if ("dataset" %in% names(cell_meta)) uniqueN(cell_meta$dataset) else NA_integer_
n_samples <- if ("cnv_sample" %in% names(cell_meta)) uniqueN(cell_meta$cnv_sample) else if ("sample_id" %in% names(cell_meta)) uniqueN(cell_meta$sample_id) else NA_integer_

perturbations <- data.table(
  perturbation = c(
    "HNF4A knockout", "HNF4A restore/OE", "PPARA knockout", "PPARA restore/OE",
    "EGR1 knockout", "CEBPB knockout", "JUN knockout", "JUND knockout", "FOS knockout",
    "AP-1 member-KO aggregate", "SOX4 knockout", "SOX4 OE",
    "HLF knockout", "IRF1 knockout", "MAFB knockout", "MAFF knockout", "MYC knockout",
    "matched random TF controls"
  ),
  tf = c(
    "HNF4A", "HNF4A", "PPARA", "PPARA", "EGR1", "CEBPB", "JUN", "JUND", "FOS",
    "AP1_AGGREGATE", "SOX4", "SOX4", "HLF", "IRF1", "MAFB", "MAFF", "MYC", "CONTROL_POOL"
  ),
  perturbation_type = c(
    "knockout", "restore_or_oe", "knockout", "restore_or_oe", rep("knockout", 5),
    "member_ko_aggregate", "knockout", "overexpression", rep("knockout", 5), "matched_controls"
  )
)
perturbations[, axis := figure6_axis_for_tf(tf)]
perturbations[tf == "AP1_AGGREGATE", axis := "stress_axis"]
perturbations[tf == "CONTROL_POOL", axis := "control"]
perturbations[, available := fifelse(
  perturbation_type == "knockout", tf %in% available_ko,
  fifelse(perturbation_type == "member_ko_aggregate", all(FIGURE6_AP1_MEMBERS %in% available_ko),
          fifelse(perturbation_type == "matched_controls", sum(FIGURE6_AXIS_TFS$control %in% available_ko) >= 1L, FALSE))
)]
perturbations[, primary_or_sensitivity := fifelse(
  tf %in% c("HNF4A", "EGR1", "SOX4") & available, "primary",
  fifelse(available, "sensitivity", "unavailable")
)]
perturbations[, signed_effect_available := available & perturbation_type %in% c("knockout", "member_ko_aggregate")]
perturbations[, magnitude_only := FALSE]
perturbations[, n_cells := fifelse(available, n_cells, NA_integer_)]
perturbations[, n_datasets := fifelse(available, n_datasets, NA_integer_)]
perturbations[, n_samples := fifelse(available, n_samples, NA_integer_)]
perturbations[, review_risk := fifelse(
  perturbation_type %in% c("restore_or_oe", "overexpression"), "not_available_no_existing_simulation",
  fifelse(perturbation_type == "member_ko_aggregate", "aggregate_of_single_gene_KOs_not_combined_perturbation",
          fifelse(perturbation_type == "matched_controls" & sum(FIGURE6_AXIS_TFS$control %in% available_ko) < 10,
                  "control_pool_below_requested_10_to_50_per_candidate", ""))
)]

checks <- rbindlist(list(
  perturbations[, .(
    check_id = "perturbation_availability", perturbation, axis,
    status = fifelse(available, "available", "unavailable"),
    primary_or_sensitivity, signed_effect_available, magnitude_only,
    n_cells, n_samples, n_datasets, review_risk,
    details = fifelse(available, "Existing CellOracle KO output or explicitly labelled aggregate", "No existing result; no value will be imputed")
  )],
  data.table(
    check_id = c(
      "delta_embedding", "predicted_delta_expression", "celloracle_target_network", "sctenifold_result",
      "sctenifold_directionality", "figure5_temporal", "external_validation", "negative_controls",
      "generic_stress_score", "proliferation_score", "hypoxia_score", "s_g2m_scores",
      "restore_oe", "ap1_combined_suppression", "figure4_legacy_mapping"
    ),
    perturbation = "global", axis = "all",
    status = c(
      if (file.exists(paths$celloracle_cells) && file.exists(paths$celloracle_grid)) "available" else "unavailable",
      if (file.exists(paths$oracle_object)) "reconstructable_from_frozen_oracle" else "partial_top50_only",
      if (file.exists(paths$celloracle_grn)) "available" else "unavailable",
      if (file.exists(paths$sct_main) && nrow(sct_contract)) "available" else "unavailable",
      "magnitude_only", if (file.exists(paths$figure5_temporal)) "available" else "unavailable",
      if (file.exists(paths$external_validation)) "available" else "unavailable",
      if (sum(FIGURE6_AXIS_TFS$control %in% available_ko) >= 1) "available_limited_pool" else "unavailable",
      if ("driver_main_strict__module_Stressed_Injured" %in% names(cell_meta)) "available" else if (file.exists(paths$oracle_object)) "reconstructable_from_frozen_oracle" else "unavailable",
      if ("proliferation_score_z" %in% names(cell_meta)) "available" else if (file.exists(paths$oracle_object)) "reconstructable_from_frozen_oracle" else "unavailable",
      "unavailable", "unavailable", "unavailable", "unavailable",
      if (file.exists(paths$three_axis_audit)) "available_legacy_figure2_filename" else "unavailable"
    ),
    primary_or_sensitivity = "audit", signed_effect_available = FALSE, magnitude_only = FALSE,
    n_cells = NA_integer_, n_samples = NA_integer_, n_datasets = NA_integer_,
    review_risk = c(
      "", "full_frozen_signature_requires_figure6_raw_export", "", "existing_module7_main_strict_low_replication",
      "distance_Z_FC_not_interpreted_as_signed_without_method_validation", "temporal_final_label_temporal_not_supported",
      "", "fewer_than_10_matched_controls_per_candidate", "generic_stress_proxy_only", "",
      "hypoxia_adjustment_not_testable", "separate_S_and_G2M_adjustment_not_testable",
      "no_restore_or_OE_outputs", "AP1_is_member_KO_aggregate_only", "SOX4 Figure 4 is stored under legacy Figure 2 names"
    ),
    details = c(
      "Cell-level saved delta_embedding and shared grids exist",
      "Fitted Oracle exists; Figure 6 exporter will compute frozen-programme delta expression without embedding",
      "Filtered CellOracle links exist", "All existing Module 7 contracts are audited",
      "Only distance/rank/overlap/pathway magnitude comparisons are permitted",
      "Module 9.1 timing outputs are readable", "Module 8 external report is readable",
      "Available simulated controls are used as a limited empirical pool",
      "Module 5 stressed/injured score is the available generic-stress proxy",
      "Existing proliferation score is available", "No frozen hypoxia score was found",
      "No separate frozen S and G2M scores were found", "KO signs will not be inverted to fabricate restore",
      "JUN/JUNB/JUND/FOS/ATF3 KOs will be summarized by a robust median",
      "Three-axis audit explicitly maps the SOX4 legacy implementation to Figure 4"
    )
  )
), fill = TRUE)

input_manifest <- data.table(
  input_type = names(paths),
  file_path = vapply(paths, figure6_norm_path, character(1))
)
input_manifest[, `:=`(
  axis = "all", perturbation = "none", method = "audit", state_or_subset = "all",
  version = "frozen_existing", source_figure = "Figure2-5", primary_or_sensitivity = "primary",
  signed_or_unsigned = "metadata", direction_interpretability = "metadata", notes = "Frozen existing project input"
)]
input_manifest[input_type %in% c("celloracle_report", "celloracle_cells", "celloracle_gene_delta", "celloracle_tf_delta", "celloracle_grid", "celloracle_grn", "celloracle_metadata", "oracle_object"),
  `:=`(perturbation = "all KO", method = "CellOracle", state_or_subset = "all CellOracle states", version = "module6_frozen", source_figure = "Figure2-4")]
input_manifest[input_type == "celloracle_cells", `:=`(signed_or_unsigned = "signed_embedding", direction_interpretability = "signed")]
input_manifest[input_type == "celloracle_gene_delta", `:=`(signed_or_unsigned = "signed_expression_partial", direction_interpretability = "signed", notes = "Top-50 state means; full frozen modules reconstructed from fitted Oracle")]
input_manifest[input_type == "celloracle_tf_delta", `:=`(signed_or_unsigned = "signed_TF_delta", direction_interpretability = "signed")]
input_manifest[input_type == "celloracle_grid", `:=`(signed_or_unsigned = "signed_embedding", direction_interpretability = "signed")]
input_manifest[input_type == "celloracle_grn", `:=`(signed_or_unsigned = "directed_GRN", direction_interpretability = "directed_network")]
input_manifest[input_type == "oracle_object", `:=`(signed_or_unsigned = "object", direction_interpretability = "raw_object", notes = "Fitted Oracle reused without changing frozen KO parameters")]
input_manifest[input_type %in% c("sct_main", "sct_contract"), `:=`(perturbation = "all KO", method = "scTenifoldKnk", state_or_subset = "main_strict", version = "module7.2_frozen", signed_or_unsigned = "unsigned", direction_interpretability = "magnitude_only", notes = "Manifold distance is unsigned")]
input_manifest[input_type == "scenic_regulons", `:=`(method = "pySCENIC", version = "module6.3", signed_or_unsigned = "regulon", direction_interpretability = "activity")]
input_manifest[input_type == "cistarget_regulons", `:=`(method = "cisTarget", version = "module6.3c", signed_or_unsigned = "regulon", direction_interpretability = "activity")]
input_manifest[input_type == "cellrank_fate", `:=`(method = "CellRank", state_or_subset = "main_strict", version = "module6.2", signed_or_unsigned = "probability", direction_interpretability = "fate_probability")]
input_manifest[input_type == "strict_cnv", `:=`(method = "CNV", state_or_subset = "strict", version = "module6.1", signed_or_unsigned = "label", direction_interpretability = "baseline_only")]
input_manifest[input_type %in% c("figure5_temporal", "figure5_order"), `:=`(method = "trajectory", version = "module9.1", source_figure = "Figure5", signed_or_unsigned = "temporal", direction_interpretability = "temporal")]
input_manifest[input_type %in% c("module9_direction", "module9_restore"), `:=`(method = "network integration", version = "module9.2", source_figure = "Figure5", signed_or_unsigned = "mixed", direction_interpretability = "guardrailed")]
input_manifest[input_type == "module9_pseudobulk", `:=`(method = "path analysis", version = "module9.3", source_figure = "Figure5", state_or_subset = "sample_pseudobulk")]
input_manifest[input_type == "external_validation", `:=`(method = "external validation", version = "module8", source_figure = "Figure5")]
input_manifest[input_type == "frozen_module_genes", `:=`(method = "trajectory signatures", version = "module5.4", source_figure = "Figure5", notes = "Frozen programme genes")]
input_manifest[input_type == "hnf4a_vector", `:=`(axis = "A", perturbation = "HNF4A KO", method = "CellOracle", state_or_subset = "identity-high", version = "Figure2", source_figure = "Figure2", signed_or_unsigned = "signed", direction_interpretability = "signed_KO", notes = "Formal HNF4A KO field")]
input_manifest[input_type == "egr1_vector", `:=`(axis = "B", perturbation = "EGR1 KO", method = "CellOracle", state_or_subset = "stress-transition", version = "Figure3", source_figure = "Figure3", signed_or_unsigned = "signed", direction_interpretability = "signed_KO", notes = "Formal EGR1 KO field")]
input_manifest[input_type == "sox4_vector", `:=`(axis = "C", perturbation = "SOX4 KO", method = "CellOracle", state_or_subset = "malignant-like", version = "Figure4_legacy", source_figure = "Figure4", signed_or_unsigned = "signed", direction_interpretability = "signed_KO", notes = "Formal SOX4 KO field under legacy Figure 2 filename")]
input_manifest[input_type == "three_axis_audit", `:=`(method = "audit", version = "three_axis_audit", source_figure = "Figure2-4", primary_or_sensitivity = "audit", notes = "Legacy mapping is explicit")]
input_manifest[, exists := file.exists(file_path)]
input_manifest[, size_bytes := fifelse(exists, file.info(file_path)$size, NA_real_)]

protected_roots <- file.path(PROJECT_ROOT, c("scripts", "metadata", "figures", "reports"))
protected <- unlist(lapply(protected_roots, function(root) {
  if (!dir.exists(root)) return(character())
  list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
}), use.names = FALSE)
protected <- protected[file.exists(protected) & !dir.exists(protected)]
rel <- substring(normalizePath(protected, mustWork = FALSE), nchar(normalizePath(PROJECT_ROOT)) + 2L)
keep <- grepl("figure[1-5]|module5|trajectory|figure2c_sox4|figure2d_sox4", rel, ignore.case = TRUE) &
  !grepl("figure6", rel, ignore.case = TRUE)
protected <- protected[keep]
rel <- rel[keep]
info <- file.info(protected)
snapshot <- data.table(
  file_path = gsub("\\\\", "/", rel),
  size_bytes = info$size,
  modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
  md5 = NA_character_
)
small <- which(is.finite(snapshot$size_bytes) & snapshot$size_bytes <= 200 * 1024^2)
if (length(small)) snapshot$md5[small] <- unname(tools::md5sum(protected[small]))

before_path <- file.path(FIGURE6_METADATA_DIR, "figure6_protected_assets_before.tsv")
if (!file.exists(before_path) || identical(figure6_get_arg("--refresh-protected-baseline", "false"), "true")) {
  figure6_fwrite(snapshot, before_path)
}

figure6_fwrite(input_manifest, file.path(FIGURE6_METADATA_DIR, "figure6_frozen_input_manifest.tsv"))
figure6_fwrite(checks, file.path(FIGURE6_METADATA_DIR, "figure6_preflight_report.tsv"))
palette_report <- figure6_palette_contract()
figure6_write_json(palette_report, file.path(FIGURE6_METADATA_DIR, "figure6_palette_report.json"))

package_table <- figure6_package_versions(c(
  "ggplot2", "ggsci", "patchwork", "cowplot", "data.table", "dplyr", "ComplexUpset",
  "igraph", "ggraph", "lavaan", "boot", "rsample", "jsonlite", "broom", "scales", "ggrepel", "tidygraph"
))
figure6_fwrite(package_table, file.path(FIGURE6_METADATA_DIR, "figure6_r_package_versions.tsv"))

preflight_report <- list(
  module = "Figure 6 preflight and frozen input audit",
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  palette = palette_report,
  n_available_knockouts = length(available_ko),
  available_knockouts = as.list(available_ko),
  unavailable_restore_or_oe = as.list(perturbations[perturbation_type %in% c("restore_or_oe", "overexpression") & !available, perturbation]),
  ap1_status = "AP-1 member-KO aggregate; not a combined knockout",
  sctenifold_directionality = "magnitude_only",
  n_preflight_rows = nrow(checks),
  n_unavailable_or_limited = checks[status %in% c("unavailable", "magnitude_only", "available_limited_pool", "partial_top50_only"), .N],
  protected_asset_count = nrow(snapshot),
  outputs = list(
    tsv = figure6_norm_path(file.path(FIGURE6_METADATA_DIR, "figure6_preflight_report.tsv")),
    frozen_manifest = figure6_norm_path(file.path(FIGURE6_METADATA_DIR, "figure6_frozen_input_manifest.tsv")),
    protected_baseline = figure6_norm_path(before_path),
    package_versions = figure6_norm_path(file.path(FIGURE6_METADATA_DIR, "figure6_r_package_versions.tsv"))
  ),
  review_risk_flags = as.list(unique(checks[nzchar(review_risk), review_risk]))
)
figure6_write_json(preflight_report, file.path(FIGURE6_METADATA_DIR, "figure6_preflight_report.json"))

cat(jsonlite::toJSON(list(
  status = "preflight_complete",
  available_knockouts = available_ko,
  restore_or_oe_available = FALSE,
  ap1_combined_suppression_available = FALSE,
  sctenifold_signed_direction_available = FALSE,
  protected_assets = nrow(snapshot),
  palette = unname(lancet_palette)
), auto_unbox = TRUE, pretty = TRUE), "\n")
