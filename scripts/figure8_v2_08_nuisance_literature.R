#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))
source(file.path("scripts", "figure8_v2_04_fetch_external_resources.R"))

figure8_v2_curated_nuisance_flag <- function(moa) {
  text <- tolower(as.character(moa %||% ""))
  as.numeric(grepl("protein synthesis|translation|ribosom|tubulin|microtubule|antimitotic|topoisomerase|dna damage|dna intercalat|mitochond|apoptosis inducer", text))
}

figure8_v2_combined_nuisance_penalty <- function(model_penalty, moa, prism_class) {
  values <- c(
    ifelse(is.finite(model_penalty), model_penalty, NA_real_),
    0.90 * figure8_v2_curated_nuisance_flag(moa),
    ifelse(!is.na(prism_class) && prism_class == "broad_cytotoxicity", 1, ifelse(!is.na(prism_class) && prism_class == "pan_cancer_activity", 0.80, 0))
  )
  if (!any(is.finite(values))) return(NA_real_)
  max(values, na.rm = TRUE)
}

figure8_v2_literature_direction <- function(text) {
  text <- tolower(as.character(text %||% ""))
  if (grepl("did not|no significant|no overall survival|no os benefit|failed|lack of efficacy|not improve|terminated|toxicity|adverse event", text)) return("negative")
  if (grepl("inhibit|suppression|antitumor|anti-tumor|xenograft|sensitive|potentiat|reduced tumor|reduced tumour", text)) return("positive_preclinical")
  "contextual"
}

figure8_v2_literature_primary_score_weight <- function() 0

figure8_v2_literature_candidate_set <- function(...) {
  values <- tolower(as.character(unlist(list(...), use.names = FALSE)))
  sort(unique(values[!is.na(values) & nzchar(values)]))
}

figure8_v2_europepmc_search <- function(candidate, page_size = 3L) {
  query <- paste0('"', candidate, '" AND (hepatocellular OR "liver cancer")')
  url <- paste0(
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
    utils::URLencode(query, reserved = TRUE), "&format=json&pageSize=", page_size
  )
  response <- tryCatch(figure8_v2_fetch_json(url), error = function(e) e)
  if (inherits(response, "error")) {
    return(data.table(candidate = candidate, evidence_type = "literature_search_audit", evidence_direction = "query_failed", title = NA_character_, abstract = NA_character_, PMID = NA_character_, PMCID = NA_character_, DOI = NA_character_, source_url = url, query_status = conditionMessage(response)))
  }
  results <- tryCatch(as.data.table(response$resultList$result), error = function(e) data.table())
  if (!nrow(results)) {
    return(data.table(candidate = candidate, evidence_type = "literature_search_audit", evidence_direction = "not_found", title = NA_character_, abstract = NA_character_, PMID = NA_character_, PMCID = NA_character_, DOI = NA_character_, source_url = url, query_status = "completed_no_results"))
  }
  for (nm in c("title", "abstractText", "pmid", "pmcid", "doi")) if (!nm %in% names(results)) results[, (nm) := NA_character_]
  results[, .(
    candidate = candidate,
    evidence_type = "Europe_PMC_exact_name_HCC_liver_search",
    evidence_direction = vapply(paste(title, abstractText), figure8_v2_literature_direction, character(1)),
    title = as.character(title), abstract = as.character(abstractText),
    PMID = as.character(pmid), PMCID = as.character(pmcid), DOI = as.character(doi),
    source_url = fifelse(!is.na(pmid), paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/"), url),
    query_status = "completed"
  )]
}

figure8_v2_manual_literature <- function() {
  data.table(
    candidate = c("everolimus", "everolimus", "everolimus", "tasquinimod", "tasquinimod", "tasquinimod", "emetine", "cephaeline", "tipifarnib-P2"),
    evidence_type = c("preclinical_HCC", "randomized_phase3_HCC", "phase1_HCC_safety", "randomized_phase2_mCRPC", "randomized_phase3_mCRPC", "single_arm_phase2_HCC", "mechanism_cytotoxicity", "mechanism_cytotoxicity", "non_HCC_preclinical_context"),
    evidence_direction = c("positive_preclinical", "negative", "negative_safety", "positive_clinical_non_HCC", "negative", "inconclusive_single_arm", "nuisance_mechanism", "nuisance_mechanism", "contextual_non_HCC"),
    disease_context = c("HCC xenograft", "advanced HCC after sorafenib", "advanced HCC", "metastatic castration-resistant prostate cancer", "metastatic castration-resistant prostate cancer", "advanced/metastatic HCC cohort", "general cancer-cell translation inhibition", "general cancer-cell translation inhibition", "head and neck squamous carcinoma"),
    study_stage = c("preclinical", "Phase III", "Phase I", "Phase II", "Phase III", "Phase II single-arm", "mechanistic", "mechanistic", "preclinical"),
    PMID = c("18466352", "25058218", "23134470", "21931019", "27298414", NA, NA, NA, NA),
    PMCID = c(NA, NA, NA, NA, NA, NA, "PMC1370947", "PMC1370947", "PMC10543974"),
    DOI = c(NA, "10.1001/jama.2014.7189", NA, "10.1200/JCO.2011.35.6295", NA, NA, NA, NA, NA),
    NCT = c(NA, "NCT01035229", "NCT00390195", NA, "NCT01234311", "NCT01743469", NA, NA, NA),
    title = c(
      "RAD001 (everolimus) inhibits tumour growth in xenograft models of human hepatocellular carcinoma",
      "Effect of everolimus on survival in advanced hepatocellular carcinoma after failure of sorafenib: EVOLVE-1",
      "Randomised clinical trial: comparison of two everolimus dosing schedules in advanced HCC",
      "Phase II randomized study of tasquinimod in metastatic castrate-resistant prostate cancer",
      "Randomized Phase III study of tasquinimod in metastatic castration-resistant prostate cancer",
      "Tasquinimod in hepatocellular, ovarian, renal and gastric cancers",
      "Eukaryotic protein synthesis inhibitors identified by comparison of cytotoxicity profiles",
      "Eukaryotic protein synthesis inhibitors identified by comparison of cytotoxicity profiles",
      "Tipifarnib plus PI3Kalpha inhibition in HNSCC models"
    ),
    summary = c(
      "Dose-dependent growth inhibition in patient-derived HCC xenografts; preclinical only.",
      "No overall-survival benefit: HR 1.05, 95% CI 0.86-1.27, P=0.68; median OS 7.6 vs 7.3 months.",
      "Dose-limiting toxicities and hepatitis-B flare risk were observed; not efficacy confirmation.",
      "Improved progression-free outcomes in a non-HCC Phase II context.",
      "Improved radiographic PFS but no overall-survival benefit; grade >=3 adverse events were more frequent.",
      "Completed single-arm study; 53 HCC patients treated; no randomized comparator, so efficacy remains inconclusive.",
      "Emetine is a reference eukaryotic protein-synthesis inhibitor identified by cytotoxicity profiling.",
      "Cephaeline is an ipecac alkaloid protein-synthesis inhibitor identified by cytotoxicity profiling.",
      "Preclinical activity was reported in HNSCC, not HCC."
    ),
    source_url = c(
      "https://pubmed.ncbi.nlm.nih.gov/18466352/", "https://pubmed.ncbi.nlm.nih.gov/25058218/", "https://pubmed.ncbi.nlm.nih.gov/23134470/",
      "https://pubmed.ncbi.nlm.nih.gov/21931019/", "https://pubmed.ncbi.nlm.nih.gov/27298414/", "https://clinicaltrials.gov/study/NCT01743469",
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC1370947/", "https://pmc.ncbi.nlm.nih.gov/articles/PMC1370947/", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543974/"
    ),
    query_status = "manually_verified_primary_or_official_source"
  )
}

figure8_v2_nuisance_literature <- function() {
  v1_penalty <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8_toxicity_stress_penalty.tsv"))
  annotation <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_compound_moa_target_annotation.tsv"))
  prism <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_prism_viability.tsv"))
  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  cross <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_concordance.tsv"))

  nuisance_cols <- c("proliferation_penalty", "generic_stress_penalty", "dna_damage_penalty", "translation_inhibition_penalty", "mitochondrial_toxicity_penalty", "pan_cytotoxicity_penalty")
  v1_penalty[, model_nuisance_penalty := apply(.SD, 1, max, na.rm = TRUE), .SDcols = nuisance_cols]
  v1_penalty[!is.finite(model_nuisance_penalty), model_nuisance_penalty := NA_real_]
  nuisance <- merge(ranking[, .(BRD_ID = compound, canonical_name, v2_primary_rank, candidate_analysis_universe)], v1_penalty, by.x = "BRD_ID", by.y = "compound", all.x = TRUE)
  nuisance <- merge(nuisance, annotation[, .(BRD_ID, curated_MoA, curated_targets, curated_annotation_available)], by = "BRD_ID", all.x = TRUE)
  nuisance <- merge(nuisance, prism[, .(BRD_ID, prism_phenotype_class, prism_gate_class)], by = "BRD_ID", all.x = TRUE)
  nuisance[, curated_nuisance_flag := vapply(curated_MoA, figure8_v2_curated_nuisance_flag, numeric(1))]
  nuisance[, nuisance_penalty := mapply(figure8_v2_combined_nuisance_penalty, model_nuisance_penalty, curated_MoA, prism_gate_class)]
  nuisance[, nuisance_interpretation := fifelse(nuisance_penalty >= 0.75, "nonspecific_apparent_reversal", fifelse(is.na(nuisance_penalty), "unavailable", "below_high_nuisance_gate"))]
  nuisance[, `:=`(
    normal_cell_safety_established = FALSE,
    penalty_source = "frozen local gene-set DrugReflector controls + curated MoA + PRISM cancer-cell phenotype",
    broad_transcriptional_suppression_status = fifelse(is.na(broad_transcriptional_suppression_status), "unknown", broad_transcriptional_suppression_status)
  )]
  figure8_v2_write_tsv(nuisance, "figure8_v2_nuisance_penalties.tsv")

  gene_sets <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8_nuisance_gene_sets.tsv"))
  gene_set_audit <- gene_sets[, .(
    n_unique_genes = uniqueN(gene), n_terms = uniqueN(term), sources = paste(sort(unique(source)), collapse = ";"),
    source_status = "frozen_local_gene_set"
  ), by = nuisance_set]
  figure8_v2_write_tsv(gene_set_audit, "figure8_v2_nuisance_gene_set_source_audit.tsv")

  top20 <- ranking[order(v2_primary_rank)][1:20, unique(tolower(canonical_name))]
  cross_names <- unique(tolower(cross[n_support_frameworks >= 2 & !is.na(canonical_name), canonical_name]))
  historical <- unique(tolower(ranking[historical_v1_top_or_reference == TRUE, canonical_name]))
  prism_names <- unique(tolower(prism[prism_gate_class == "hcc_liver_enriched", canonical_name]))
  integrated_path <- file.path(FIGURE8_V2_METADATA, "figure8_v2_integrated_candidate_evidence.tsv")
  integrated_names <- if (file.exists(integrated_path)) {
    integrated <- fread(integrated_path)
    integrated[candidate_analysis_universe == TRUE][order(candidate_priority_rank)][1:20, canonical_name]
  } else character()
  candidates <- figure8_v2_literature_candidate_set(top20, cross_names, historical, prism_names, integrated_names)
  automatic <- rbindlist(lapply(candidates, figure8_v2_europepmc_search), fill = TRUE)
  manual <- figure8_v2_manual_literature()
  literature <- rbindlist(list(manual, automatic), fill = TRUE)
  literature[, `:=`(
    retrieval_date = "2026-08-18",
    primary_computational_score_weight = figure8_v2_literature_primary_score_weight(),
    literature_role = "contextual_evidence_only"
  )]
  setorder(literature, candidate, evidence_direction, evidence_type)
  figure8_v2_write_tsv(literature, "figure8_v2_candidate_literature_evidence.tsv")

  figure8_v2_write_json(list(
    module = "figure8_v2_nuisance_literature", status = "completed",
    high_nuisance_candidate_count = nuisance[nuisance_penalty >= 0.75, .N],
    candidate_literature_queries = length(candidates),
    literature_rows = nrow(literature),
    literature_primary_score_weight = 0,
    safety_boundary = "Nuisance scores and PRISM do not establish normal-cell safety."
  ), "figure8_v2_nuisance_literature_report.json")
  invisible(list(nuisance = nuisance, literature = literature))
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_nuisance_literature()
  cat("FIGURE8_V2_NUISANCE high=", sum(result$nuisance$nuisance_penalty >= 0.75, na.rm = TRUE), " literature_rows=", nrow(result$literature), "\n", sep = "")
}
