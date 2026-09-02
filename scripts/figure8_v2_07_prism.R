#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))
source(file.path("scripts", "figure8_v2_06_moa_network.R"))
source(file.path("scripts", "figure8_v2_04_fetch_external_resources.R"))

figure8_v2_context_label <- function(cell_line_name, oncotree_code, primary_disease, lineage) {
  name <- toupper(gsub("[^A-Z0-9]", "", as.character(cell_line_name)))
  code <- toupper(as.character(oncotree_code))
  disease <- tolower(as.character(primary_disease))
  lineage_lower <- tolower(as.character(lineage))
  if (name == "HCC515") return("lung_adenocarcinoma_non_liver")
  if (name == "HA1E") return("kidney_derived_non_liver")
  if (name == "HEPG2") return("liver_derived_hepatoblastoma_like_caveat")
  if (!is.na(code) && code == "HCC") return("adult_hepatocellular_carcinoma")
  if ((!is.na(code) && code == "LIHB") || grepl("hepatoblastoma", disease)) return("liver_derived_hepatoblastoma_like_caveat")
  if ((!is.na(code) && code == "HCCIHCH") || grepl("hepatocellular carcinoma plus", disease)) return("mixed_hcc_intrahepatic_cholangiocarcinoma")
  if (!is.na(lineage_lower) && lineage_lower == "liver") return("other_liver_derived_cancer")
  if (grepl("non-cancerous", disease)) return("non_cancerous_non_liver")
  "other_cancer_non_liver"
}

figure8_v2_choose_secondary_screen <- function(screen_ids) {
  screen_ids <- unique(as.character(screen_ids[!is.na(screen_ids)]))
  if ("MTS010" %in% screen_ids) "MTS010" else if ("HTS002" %in% screen_ids) "HTS002" else if (length(screen_ids)) sort(screen_ids)[[1]] else NA_character_
}

figure8_v2_add_secondary_sensitivity <- function(x) {
  x <- as.data.table(copy(x))
  x[, secondary_sensitivity := -as.numeric(auc)]
  x
}

figure8_v2_combine_prism_classes <- function(primary_class, secondary_class) {
  classes <- c(primary_class, secondary_class)
  classes <- classes[!is.na(classes) & !classes %in% c("unavailable", "insufficient_liver_lines")]
  if (!length(classes)) return("unavailable")
  has_enriched <- any(classes == "hcc_liver_enriched")
  has_broad <- any(classes == "broad_cytotoxicity")
  if (has_enriched && has_broad) return("discordant_phenotype")
  if (sum(classes == "hcc_liver_enriched") >= 2L) return("hcc_liver_enriched_replicated")
  if (has_enriched) return("hcc_liver_enriched_single_release")
  if (has_broad) return("broad_cytotoxicity")
  if (any(classes == "pan_cancer_activity")) return("pan_cancer_activity")
  "no_enriched_support"
}

figure8_v2_cellosaurus_disease <- function(record) {
  tryCatch({
    cell <- record$Cellosaurus$`cell-line-list`[[1]]
    disease <- cell$`disease-list`
    if (is.data.frame(disease) && "label" %in% names(disease)) paste(unique(disease$label), collapse = ";") else NA_character_
  }, error = function(e) NA_character_)
}

figure8_v2_prism <- function() {
  raw_dir <- file.path(FIGURE8_V2_DATA, "external_raw")
  model <- fread(file.path(raw_dir, "Model.csv"))
  model[, cell_line := fifelse(!is.na(StrippedCellLineName) & StrippedCellLineName != "", StrippedCellLineName, gsub("[^A-Za-z0-9]", "", CellLineName))]
  model[, verified_context := mapply(figure8_v2_context_label, CellLineName, OncotreeCode, OncotreePrimaryDisease, OncotreeLineage)]
  model[, context_source := "DepMap 23Q2 Public v4 Model.csv"]

  prism_cells <- fread(file.path(raw_dir, "Repurposing_Public_23Q2_Cell_Line_Meta_Data.csv"))
  primary_depmap <- unique(prism_cells$depmap_id)
  secondary_cells <- fread(file.path(raw_dir, "secondary-screen-cell-line-info.csv"))
  secondary_depmap <- unique(secondary_cells$depmap_id)
  key_ids <- c("ACH-000872", "ACH-001310", "ACH-000739")
  cell_audit <- model[ModelID %in% unique(c(primary_depmap, secondary_depmap, key_ids)), .(
    depmap_id = ModelID, cell_line, CellLineName, CCLEName, RRID, SourceType,
    sample_collection_site = SampleCollectionSite, oncotree_code = OncotreeCode,
    oncotree_subtype = OncotreeSubtype, primary_disease = OncotreePrimaryDisease,
    lineage = OncotreeLineage, verified_context, context_source,
    in_prism_23q2 = ModelID %in% primary_depmap, in_prism_19q4_secondary = ModelID %in% secondary_depmap
  )]

  cellosaurus_dir <- file.path(raw_dir, "cellosaurus_key_records")
  dir.create(cellosaurus_dir, recursive = TRUE, showWarnings = FALSE)
  key_records <- data.table(
    depmap_id = key_ids,
    cell_line = c("HCC515", "HA1E", "HEPG2"),
    cellosaurus_accession = c("CVCL_5136", "CVCL_VU89", "CVCL_0027")
  )
  key_records[, cellosaurus_url := paste0("https://api.cellosaurus.org/cell-line/", cellosaurus_accession, "?format=json")]
  key_records[, cellosaurus_disease := NA_character_]
  key_records[, cellosaurus_status := "not_retrieved"]
  for (idx in seq_len(nrow(key_records))) {
    record <- tryCatch(figure8_v2_fetch_json(key_records$cellosaurus_url[idx]), error = function(e) e)
    if (!inherits(record, "error")) {
      write_json(record, file.path(cellosaurus_dir, paste0(key_records$cellosaurus_accession[idx], ".json")), pretty = TRUE, auto_unbox = TRUE)
      key_records[idx, `:=`(cellosaurus_disease = figure8_v2_cellosaurus_disease(record), cellosaurus_status = "verified")]
    } else key_records[idx, cellosaurus_status := paste0("failed: ", conditionMessage(record))]
  }
  cell_audit <- merge(cell_audit, key_records[, .(depmap_id, cellosaurus_accession, cellosaurus_url, cellosaurus_disease, cellosaurus_status)], by = "depmap_id", all.x = TRUE)
  figure8_v2_write_tsv(cell_audit, "figure8_v2_cell_line_metadata_audit.tsv")

  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  ranking[, brd_core := figure8_v2_brd_core(compound)]
  candidate_cores <- unique(ranking[candidate_analysis_universe == TRUE, brd_core])

  primary <- fread(file.path(raw_dir, "Repurposing_Public_23Q2_LFC_COLLAPSED.csv"))
  primary[, `:=`(depmap_id = sub("::.*$", "", row_id), brd_core = figure8_v2_brd_core(broad_id), primary_sensitivity = -as.numeric(LFC))]
  primary <- primary[brd_core %in% candidate_cores]
  primary <- merge(primary, cell_audit[, .(depmap_id, cell_line, verified_context, primary_disease, lineage)], by = "depmap_id", all.x = TRUE)
  primary[, cancer_model := !verified_context %in% c("non_cancerous_non_liver") & !is.na(verified_context)]
  figure8_v2_write_tsv(primary, "figure8_v2_prism_primary_candidate_cell_values.tsv.gz", compress = TRUE)

  primary_summary <- primary[cancer_model == TRUE, {
    adult <- primary_sensitivity[verified_context == "adult_hepatocellular_carcinoma"]
    hepatoblastoma <- primary_sensitivity[verified_context == "liver_derived_hepatoblastoma_like_caveat"]
    mixed <- primary_sensitivity[verified_context == "mixed_hcc_intrahepatic_cholangiocarcinoma"]
    liver <- primary_sensitivity[verified_context %in% c("adult_hepatocellular_carcinoma", "liver_derived_hepatoblastoma_like_caveat", "mixed_hcc_intrahepatic_cholangiocarcinoma", "other_liver_derived_cancer")]
    non_liver <- primary_sensitivity[!verified_context %in% c("adult_hepatocellular_carcinoma", "liver_derived_hepatoblastoma_like_caveat", "mixed_hcc_intrahepatic_cholangiocarcinoma", "other_liver_derived_cancer")]
    adult_median <- if (length(adult)) median(adult, na.rm = TRUE) else NA_real_
    .(
      primary_pan_cancer_median_sensitivity = median(primary_sensitivity, na.rm = TRUE),
      primary_adult_hcc_median_sensitivity = adult_median,
      primary_hepatoblastoma_median_sensitivity = if (length(hepatoblastoma)) median(hepatoblastoma, na.rm = TRUE) else NA_real_,
      primary_mixed_hcc_ihcc_median_sensitivity = if (length(mixed)) median(mixed, na.rm = TRUE) else NA_real_,
      primary_all_liver_median_sensitivity = if (length(liver)) median(liver, na.rm = TRUE) else NA_real_,
      primary_non_liver_median_sensitivity = if (length(non_liver)) median(non_liver, na.rm = TRUE) else NA_real_,
      primary_adult_hcc_within_cellline_percentile = if (is.finite(adult_median)) mean(primary_sensitivity <= adult_median, na.rm = TRUE) else NA_real_,
      primary_liver_minus_non_liver = if (length(adult) && length(non_liver)) median(adult, na.rm = TRUE) - median(non_liver, na.rm = TRUE) else NA_real_,
      n_primary_total_cancer_lines = uniqueN(depmap_id),
      n_primary_adult_hcc_lines = uniqueN(depmap_id[verified_context == "adult_hepatocellular_carcinoma"]),
      n_primary_hepatoblastoma_lines = uniqueN(depmap_id[verified_context == "liver_derived_hepatoblastoma_like_caveat"]),
      n_primary_all_liver_lines = uniqueN(depmap_id[verified_context %in% c("adult_hepatocellular_carcinoma", "liver_derived_hepatoblastoma_like_caveat", "mixed_hcc_intrahepatic_cholangiocarcinoma", "other_liver_derived_cancer")])
    )
  }, by = brd_core]
  primary_summary[, primary_pan_cancer_percentile := frank(primary_pan_cancer_median_sensitivity, ties.method = "average") / .N]
  primary_summary[, prism_primary_class := mapply(
    figure8_v2_prism_class,
    primary_adult_hcc_within_cellline_percentile, primary_pan_cancer_percentile,
    primary_liver_minus_non_liver, n_primary_adult_hcc_lines,
    MoreArgs = list(enrichment_delta = 0.10)
  )]
  primary_summary[, prism_primary_score := mapply(
    figure8_v2_prism_score,
    primary_adult_hcc_within_cellline_percentile, primary_pan_cancer_percentile,
    primary_liver_minus_non_liver
  )]

  secondary <- fread(file.path(raw_dir, "secondary-screen-dose-response-curve-parameters.csv"))
  secondary[, brd_core := figure8_v2_brd_core(broad_id)]
  secondary <- secondary[brd_core %in% candidate_cores & passed_str_profiling == TRUE & is.finite(as.numeric(auc)) & as.numeric(r2) >= 0.70]
  secondary[, chosen_screen := figure8_v2_choose_secondary_screen(screen_id), by = .(brd_core, depmap_id)]
  secondary <- secondary[screen_id == chosen_screen]
  secondary <- merge(secondary, cell_audit[, .(depmap_id, cell_line, verified_context, primary_disease, lineage)], by = "depmap_id", all.x = TRUE)
  secondary <- figure8_v2_add_secondary_sensitivity(secondary)
  figure8_v2_write_tsv(secondary, "figure8_v2_prism_secondary_candidate_curve_values.tsv.gz", compress = TRUE)

  secondary_summary <- secondary[, {
    adult_sensitivity <- secondary_sensitivity[verified_context == "adult_hepatocellular_carcinoma"]
    non_liver_sensitivity <- secondary_sensitivity[!verified_context %in% c("adult_hepatocellular_carcinoma", "liver_derived_hepatoblastoma_like_caveat", "mixed_hcc_intrahepatic_cholangiocarcinoma", "other_liver_derived_cancer")]
    adult_median_sensitivity <- if (length(adult_sensitivity)) median(adult_sensitivity, na.rm = TRUE) else NA_real_
    .(
    secondary_pan_cancer_median_auc = median(as.numeric(auc), na.rm = TRUE),
    secondary_pan_cancer_median_sensitivity = median(secondary_sensitivity, na.rm = TRUE),
    secondary_adult_hcc_median_auc = if (length(adult_sensitivity)) median(as.numeric(auc)[verified_context == "adult_hepatocellular_carcinoma"], na.rm = TRUE) else NA_real_,
    secondary_adult_hcc_median_sensitivity = adult_median_sensitivity,
    secondary_hepatoblastoma_median_auc = if (any(verified_context == "liver_derived_hepatoblastoma_like_caveat")) median(as.numeric(auc)[verified_context == "liver_derived_hepatoblastoma_like_caveat"], na.rm = TRUE) else NA_real_,
    secondary_adult_hcc_within_cellline_percentile = if (is.finite(adult_median_sensitivity)) mean(secondary_sensitivity <= adult_median_sensitivity, na.rm = TRUE) else NA_real_,
    secondary_liver_minus_non_liver = if (length(adult_sensitivity) && length(non_liver_sensitivity)) median(adult_sensitivity, na.rm = TRUE) - median(non_liver_sensitivity, na.rm = TRUE) else NA_real_,
    secondary_valid_ic50_median = if (any(is.finite(as.numeric(ic50)))) median(as.numeric(ic50)[is.finite(as.numeric(ic50))], na.rm = TRUE) else NA_real_,
    secondary_valid_ic50_count = sum(is.finite(as.numeric(ic50))),
    secondary_curve_lower_limit_median = median(as.numeric(lower_limit), na.rm = TRUE),
    secondary_median_r2 = median(as.numeric(r2), na.rm = TRUE),
    n_secondary_total_lines = uniqueN(depmap_id),
    n_secondary_adult_hcc_lines = uniqueN(depmap_id[verified_context == "adult_hepatocellular_carcinoma"]),
    secondary_screen_ids = paste(sort(unique(screen_id)), collapse = ";")
  )}, by = brd_core]
  secondary_summary[, secondary_pan_cancer_percentile := frank(secondary_pan_cancer_median_sensitivity, ties.method = "average") / .N]
  secondary_summary[, prism_secondary_class := mapply(
    figure8_v2_prism_class,
    secondary_adult_hcc_within_cellline_percentile, secondary_pan_cancer_percentile,
    secondary_liver_minus_non_liver, n_secondary_adult_hcc_lines,
    MoreArgs = list(enrichment_delta = 0.10)
  )]
  secondary_summary[, prism_secondary_score := mapply(
    figure8_v2_prism_score,
    secondary_adult_hcc_within_cellline_percentile, secondary_pan_cancer_percentile,
    secondary_liver_minus_non_liver
  )]

  result <- merge(ranking[, .(BRD_ID = compound, brd_core, canonical_name, inchi_key, v2_primary_rank, candidate_analysis_universe)], primary_summary, by = "brd_core", all.x = TRUE)
  result <- merge(result, secondary_summary, by = "brd_core", all.x = TRUE)
  result[, `:=`(
    primary_screen_available = !is.na(primary_pan_cancer_median_sensitivity),
    secondary_dose_response_available = !is.na(secondary_pan_cancer_median_auc),
    prism_status = fifelse(!is.na(primary_pan_cancer_median_sensitivity) & !is.na(secondary_pan_cancer_median_auc), "primary_and_secondary", fifelse(!is.na(primary_pan_cancer_median_sensitivity), "primary_only", fifelse(!is.na(secondary_pan_cancer_median_auc), "secondary_only", "not_tested"))),
    normal_cell_safety_established = FALSE,
    safety_statement = "PRISM is a cancer-cell viability assay; normal-cell safety is not established",
    primary_source = "PRISM Repurposing Public 23Q2 v4, 2.5 uM, 5-day single-dose LFC",
    secondary_source = "PRISM Repurposing 19Q4 v4 secondary dose-response; valid curve R2>=0.70; MTS010 preferred"
  )]
  result[, prism_phenotype_class := mapply(figure8_v2_combine_prism_classes, prism_primary_class, prism_secondary_class)]
  result[, prism_gate_class := fifelse(grepl("^hcc_liver_enriched", prism_phenotype_class), "hcc_liver_enriched", fifelse(prism_phenotype_class == "broad_cytotoxicity", "broad_cytotoxicity", prism_phenotype_class))]
  result[, prism_phenotype_score := rowMeans(cbind(prism_primary_score, prism_secondary_score), na.rm = TRUE)]
  result[!is.finite(prism_phenotype_score), prism_phenotype_score := NA_real_]
  figure8_v2_write_tsv(result, "figure8_v2_prism_viability.tsv")
  availability <- result[, .N, by = .(candidate_analysis_universe, prism_status)][order(-candidate_analysis_universe, prism_status)]
  availability[, `:=`(
    resource = "PRISM",
    status_definition = "not_tested remains NA and is not negative evidence",
    normal_cell_safety_established = FALSE
  )]
  figure8_v2_write_tsv(availability, "figure8_v2_prism_availability_report.tsv")
  figure8_v2_write_json(list(
    module = "figure8_v2_prism", status = "completed",
    primary_release = "Repurposing Public 23Q2 v4", secondary_release = "PRISM Repurposing 19Q4 v4 secondary",
    candidate_primary_coverage = result[candidate_analysis_universe == TRUE, mean(primary_screen_available)],
    candidate_secondary_coverage = result[candidate_analysis_universe == TRUE, mean(secondary_dose_response_available)],
    hcc_liver_enriched_count = result[prism_gate_class == "hcc_liver_enriched", .N],
    broad_cytotoxicity_count = result[prism_gate_class == "broad_cytotoxicity", .N],
    context_corrections = list(HCC515 = "lung adenocarcinoma / non-liver", HA1E = "kidney-derived / non-liver", HepG2 = "liver-derived hepatoblastoma-like caveat"),
    interpretation_boundary = "Cancer-cell viability only; no normal-cell safety or clinical efficacy inference."
  ), "figure8_v2_prism_report.json")
  invisible(result)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_prism()
  cat("FIGURE8_V2_PRISM primary=", sum(result$primary_screen_available), " secondary=", sum(result$secondary_dose_response_available), " hcc_enriched=", sum(result$prism_gate_class == "hcc_liver_enriched", na.rm = TRUE), " broad=", sum(result$prism_gate_class == "broad_cytotoxicity", na.rm = TRUE), "\n", sep = "")
}
