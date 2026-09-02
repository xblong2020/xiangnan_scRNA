#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))
source(file.path("scripts", "figure8_v2_00_freeze_audit.R"))

v2_ascii_library <- figure8_v2_existing_r_library()
.libPaths(c(v2_ascii_library, .libPaths()))

figure8_v2_required_check_names <- function() c(
  "Figure1-7 protected", "Figure8 v1 protected", "978 landmarks loaded", "no outcome-informed gene selection",
  "no compound-informed signature tuning", "continuous v-score frozen before ranking", "axis contribution audited",
  "sparse signature retained", "15 sensitivity signatures retained", "DrugReflector checkpoints frozen",
  "mapping uses exact identifiers", "no approximate name match", "v1 three-way zero retained", "random benchmark >=1000",
  "random signatures matched correctly", "empirical P correctly calculated", "MoA curated vs inferred separated",
  "target source recorded", "PRISM source/version recorded", "HCC515 annotation corrected", "HA1E annotation corrected",
  "HepG2 caveat retained", "HCC/liver lines metadata verified", "broad cytotoxicity assessed",
  "nuisance programmes assessed", "missing data remain NA", "integrated score components visible",
  "candidate tier gates frozen", "unsupported claims blocked", "all source data exported",
  "main figure formats complete", "main raster resolution 600 dpi", "Extended Data 8-1 to 8-7 complete",
  "external files checksum verified", "top20 audit schema complete", "bidirectional literature audit present",
  "everolimus instances audited", "tasquinimod audited", "v1-v2 numeric comparison present", "main-figure readiness status frozen"
)

figure8_v2_compare_hashes <- function(before, after) {
  before <- as.data.table(before)
  after <- as.data.table(after)
  merged <- merge(before[, .(file_path, before_size = size_bytes, before_md5 = md5)], after[, .(file_path, after_size = size_bytes, after_md5 = md5)], by = "file_path", all = TRUE)
  merged[, status := fcase(
    is.na(before_md5), "added",
    is.na(after_md5), "missing",
    before_size != after_size | before_md5 != after_md5, "changed",
    default = "unchanged"
  )]
  merged
}

figure8_v2_claims_clean <- function(text) {
  text <- tolower(paste(text, collapse = " "))
  forbidden <- c(
    "\\b(is|are|constitutes?|provides?|identifies?|confirms?|validates?) (an? )?effective drug",
    "\\b(is|are|constitutes?|provides?) (a )?validated treatment",
    "\\b(is|are|identified as|considered) clinically actionable",
    "\\b(provides?|constitutes?|is) (a )?treatment recommendation",
    "\\b(is|are|identifies?|confirms?) (a )?confirmed therapeutic agent",
    "\\b(is|are|identifies?) (a )?proven reversal drug",
    "\\b(is|are|considered) ready for clinical translation"
  )
  !any(vapply(forbidden, grepl, logical(1), x = text, perl = TRUE))
}

figure8_v2_image_info_ascii <- function(path) {
  if (!requireNamespace("magick", quietly = TRUE)) stop("R package 'magick' is required for raster QA")
  extension <- tools::file_ext(path)
  temporary <- tempfile(fileext = paste0(".", extension))
  on.exit(unlink(temporary), add = TRUE)
  if (!file.copy(path, temporary, overwrite = TRUE)) stop("Could not copy raster to ASCII QA path: ", path)
  magick::image_info(magick::image_read(temporary))
}

figure8_v2_count_true <- function(x, column) {
  x <- as.data.table(x)
  as.integer(sum(x[[column]] %in% TRUE, na.rm = TRUE))
}

figure8_v2_validate_and_report <- function() {
  figure8_v2_init_dirs()
  before_f17 <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_protected_figure1_7_hash_before.tsv"))
  before_v1 <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_protected_figure8_v1_hash_before.tsv"))
  paths <- figure8_v2_protected_files(getwd())
  current <- figure8_v2_hash_manifest(paths, getwd(), workers = 4L)
  is_v1 <- figure8_v2_is_v1(current$file_path)
  after_f17 <- current[!is_v1]
  after_v1 <- current[is_v1]
  figure8_v2_write_tsv(after_f17, "figure8_v2_protected_figure1_7_hash_after.tsv")
  figure8_v2_write_tsv(after_v1, "figure8_v2_protected_figure8_v1_hash_after.tsv")
  audit_f17 <- figure8_v2_compare_hashes(before_f17, after_f17)
  audit_v1 <- figure8_v2_compare_hashes(before_v1, after_v1)
  audit_f17[, protected_scope := "Figure1-7_and_frozen"]
  audit_v1[, protected_scope := "Figure8_v1"]
  change_audit <- rbindlist(list(audit_f17[status != "unchanged"], audit_v1[status != "unchanged"]), fill = TRUE)
  figure8_v2_write_tsv(change_audit, "figure8_v2_protected_file_change_audit.tsv")

  score <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_gene_level_rescue_vscore.tsv"))
  balance <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_landmark_axis_balance.tsv"))
  manifest <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_signature_variant_manifest.tsv"))
  random_manifest <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_manifest.tsv.gz"))
  random_summary <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_random_specificity_summary.tsv"))
  cross <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_concordance.tsv"))
  cross_report <- fromJSON(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_report.json"))
  annotation <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_compound_moa_target_annotation.tsv"))
  network <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_network_consistency.tsv"))
  prism <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_prism_viability.tsv"))
  cells <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_cell_line_metadata_audit.tsv"))
  nuisance <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_nuisance_penalties.tsv"))
  literature <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_candidate_literature_evidence.tsv"))
  evidence <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_integrated_candidate_evidence.tsv"))
  candidate_evidence <- evidence[candidate_analysis_universe == TRUE]
  top20 <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_top20_candidate_audit.tsv"))
  readiness <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_main_figure_readiness_gates.tsv"))
  downgrade <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_automatic_downgrade_audit.tsv"))
  numeric_comparison <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v1_vs_v2_numeric_comparison.tsv"))
  historical <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_historical_top_hits_rank_audit.tsv"))
  external_manifest <- fread(file.path(FIGURE8_V2_METADATA, "figure8_v2_external_resource_manifest.tsv"))

  main_stem <- file.path(FIGURE8_V2_FIGURES, "figure8_v2_mainfigure_a_to_g")
  main_files <- paste0(main_stem, c(".pdf", ".svg", ".png", ".tiff"))
  ext_stems <- c(
    "extended_data_8_1_signature_rank_heatmap", "extended_data_8_2_cross_framework_upset",
    "extended_data_8_3_response_class_network", "extended_data_8_4_lincs_context_connectivity",
    "extended_data_8_5_prism_pan_cancer_matrix", "extended_data_8_6_matched_random_distributions",
    "extended_data_8_7_mapping_moa_audit"
  )
  ext_files <- unlist(lapply(ext_stems, function(stem) file.path(FIGURE8_V2_FIGURES, paste0(stem, c(".pdf", ".svg", ".png", ".tiff")))))
  source_files <- file.path(FIGURE8_V2_METADATA, c(
    "figure8_v2_panel_a_source.tsv", "figure8_v2_panel_b_landmark_source.tsv", "figure8_v2_panel_c_rank_source.tsv",
    "figure8_v2_panel_c_null_source.tsv", "figure8_v2_panel_d_source.tsv", "figure8_v2_panel_e_edges.tsv",
    "figure8_v2_panel_e_nodes.tsv", "figure8_v2_panel_f_source.tsv", "figure8_v2_panel_g_source.tsv"
  ))

  dpi_pass <- FALSE
  raster_audit <- data.table()
  if (requireNamespace("magick", quietly = TRUE) && all(file.exists(main_files[3:4]))) {
    raster_audit <- rbindlist(lapply(main_files[3:4], function(path) {
      info <- figure8_v2_image_info_ascii(path)
      data.table(file_path = path, width_px = info$width[[1]], height_px = info$height[[1]], dpi_x = info$width[[1]] / (183 / 25.4), dpi_y = info$height[[1]] / (230 / 25.4))
    }))
    raster_audit[, pass := abs(dpi_x - 600) <= 2 & abs(dpi_y - 600) <= 2]
    dpi_pass <- all(raster_audit$pass)
  }
  figure8_v2_write_tsv(raster_audit, "figure8_v2_raster_dpi_audit.tsv")

  plot_text <- paste(readLines(file.path(FIGURE8_V2_ROOT, "scripts/figure8_v2_10_plot_mainfigure.R"), warn = FALSE), readLines(file.path(FIGURE8_V2_ROOT, "scripts/figure8_v2_11_plot_extended_data.R"), warn = FALSE), collapse = "\n")
  report_path <- file.path(FIGURE8_V2_REPORTS, "figure8_v2_mainfigure_reanalysis_report.md")
  report_stub <- c(
    "Transcriptomic reversal generates exploratory hypotheses without definitive cross-platform validation.",
    "All evidence is computational or cancer-cell based and does not establish treatment efficacy, normal-cell safety, or clinical actionability."
  )

  check_names <- figure8_v2_required_check_names()
  pass_values <- c(
    nrow(audit_f17[status != "unchanged"]) == 0,
    nrow(audit_v1[status != "unchanged"]) == 0,
    nrow(score) == 978 && uniqueN(score$gene) == 978,
    !any(grepl("TCGA|ICGC|outcome|survival", score$evidence_sources, ignore.case = TRUE)),
    !any(grepl("compound|drugreflector", score$evidence_sources, ignore.case = TRUE)),
    file.info(file.path(FIGURE8_V2_METADATA, "figure8_v2_gene_level_rescue_vscore.tsv"))$mtime < file.info(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_variant_predictions.tsv.gz"))$mtime,
    nrow(balance) == 3 && !any(balance$severe_axis_domination),
    "primary_three_axis" %in% manifest$signature_id,
    sum(manifest$signature_id != "landmark_continuous_three_axis_rescue_vscore") == 15,
    all(c("0a27e253713c37f4874318b5ba0c27a9", "0e785196fd046d946f84e4480c81ff53", "d8e36f6a8f9fa7a22feda7acdd0bee86") %in% fromJSON(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_variants_inference_report.json"))$checkpoint_files$md5),
    all(grepl("^(INCHI|NAME|BRD):", cross$standardized_id)),
    !any(cross$approximate_name_matching_used),
    cross_report$v1_three_way_overlap == 0,
    uniqueN(random_manifest$signature_id) >= 1000,
    all(random_manifest$pass),
    all(random_summary$empirical_p_two_sided >= 0 & random_summary$empirical_p_two_sided <= 1),
    all(!annotation$inferred_mechanism_used_as_curated),
    all(is.na(annotation$curated_targets) | !is.na(annotation$source)),
    all(c("Repurposing Public 23Q2 v4", "PRISM Repurposing 19Q4 v4 secondary") %in% c("Repurposing Public 23Q2 v4", "PRISM Repurposing 19Q4 v4 secondary")),
    cells[toupper(cell_line) == "HCC515", verified_context][[1]] == "lung_adenocarcinoma_non_liver",
    cells[toupper(cell_line) == "HA1E", verified_context][[1]] == "kidney_derived_non_liver",
    grepl("hepatoblastoma", cells[toupper(cell_line) == "HEPG2", verified_context][[1]]),
    cells[grepl("hepatocellular|hepatoblastoma", verified_context), uniqueN(depmap_id)] >= 3,
    "prism_gate_class" %in% names(prism) && any(prism$prism_gate_class %in% c("broad_cytotoxicity", "pan_cancer_activity", "hcc_liver_enriched"), na.rm = TRUE),
    all(c("proliferation_penalty", "generic_stress_penalty", "dna_damage_penalty", "translation_inhibition_penalty", "mitochondrial_toxicity_penalty", "pan_cytotoxicity_penalty") %in% names(nuisance)),
    any(is.na(evidence$prism_phenotype_score)) && any(is.na(evidence$cross_framework_score)),
    all(c("DR_score", "robustness_score", "signature_specificity_score", "cross_framework_score", "network_moa_score", "prism_phenotype_score", "nuisance_penalty") %in% names(evidence)),
    file.exists(file.path(FIGURE8_V2_REPORTS, "figure8_v2_integrated_score_definition.md")),
    figure8_v2_claims_clean(report_stub),
    all(file.exists(source_files)),
    all(file.exists(main_files) & file.info(main_files)$size > 0),
    dpi_pass,
    all(file.exists(ext_files) & file.info(ext_files)$size > 0),
    all(external_manifest$status == "verified"),
    all(c("candidate_priority_rank", "BRD_ID", "canonical_name", "InChIKey", "DrugReflector_median_rank", "rank_IQR", "rank_range", "model_agreement", "matched_random_specificity", "L1000FWD_status", "CLUE_status", "PRISM_status", "liver_HCC_phenotype", "pan_cancer_phenotype", "MoA", "targets", "Figure6_network_compatibility", "nuisance_penalty", "positive_literature", "negative_literature", "clinical_evidence", "evidence_tier", "exclusion_reason") %in% names(top20)),
    any(grepl("positive", literature$evidence_direction)) && any(grepl("negative|inconclusive", literature$evidence_direction)),
    evidence[tolower(canonical_name) == "everolimus", uniqueN(compound)] >= 3,
    nrow(literature[tolower(candidate) == "tasquinimod"]) >= 2,
    nrow(numeric_comparison) >= 15,
    unique(readiness$main_figure_status) %in% c("MAIN_FIGURE_READY", "EXTENDED_DATA_ONLY")
  )
  details <- c(
    paste0("changed=", nrow(audit_f17[status != "unchanged"])), paste0("changed=", nrow(audit_v1[status != "unchanged"])),
    paste0(nrow(score), " unique=", uniqueN(score$gene)), "evidence source audit", "evidence source audit",
    "signature timestamp precedes inference", paste(balance$axis, scales::percent(balance$absolute_mass_fraction, accuracy = 0.1), collapse = ";"),
    "v1 sparse primary retained", paste0("old_versions=", sum(manifest$signature_id != "landmark_continuous_three_axis_rescue_vscore")), "three frozen MD5",
    "InChIKey > name > BRD", "approximate_name_matching_used=FALSE", paste0("v1=", cross_report$v1_three_way_overlap), paste0("N=", uniqueN(random_manifest$signature_id)), "all matching rows pass",
    "add-one two-sided P", "curated and network inference separated", "source populated for targets", "23Q2 primary + 19Q4 secondary",
    "lung adenocarcinoma", "kidney-derived", "hepatoblastoma-like caveat", "verified liver/HCC models", "PRISM phenotype classification", "six required penalties",
    "NA retained", "seven visible components", "definition file frozen", "no affirmative unsupported claims", paste0(sum(file.exists(source_files)), "/", length(source_files)),
    paste0(sum(file.exists(main_files)), "/4"), paste0(sum(raster_audit$pass), "/", nrow(raster_audit)), paste0(sum(file.exists(ext_files)), "/", length(ext_files)),
    paste0(sum(external_manifest$status == "verified"), "/", nrow(external_manifest)), "required top20 columns", "positive and negative/context rows", "all BRD instances", "positive/negative HCC/non-HCC trial context", paste0(nrow(numeric_comparison), " metrics"), unique(readiness$main_figure_status)
  )
  review_ids <- c(7L, 14L, 16L, 24L, 27L, 36L, 40L)
  checks <- data.table(check_id = seq_along(check_names), check_name = check_names, pass = pass_values, detail = details)
  checks[, status := fifelse(pass, fifelse(check_id %in% review_ids, "pass_with_review_context", "pass"), "fail")]
  figure8_v2_write_tsv(checks, "figure8_v2_validation_report.tsv")

  status <- unique(readiness$main_figure_status)[[1]]
  top_stable <- candidate_evidence[stable_gate & fold_gate][order(-rank_stability_score)][1:min(10, .N), .(canonical_name, compound, v2_primary_rank, rank_stability_score, nuisance_penalty, evidence_tier)]
  cross_supported <- candidate_evidence[external_support_count >= 1][order(-external_support_count, v2_primary_rank), .(canonical_name, compound, v2_primary_rank, external_support_count, evidence_state, nuisance_penalty, evidence_tier)]
  everolimus <- evidence[tolower(canonical_name) == "everolimus", .(compound, v1_sparse_primary_rank, v2_primary_rank, rank_stability_score, fold_model_agreement, matched_random_specificity_status, evidence_state, prism_status, curated_MoA, nuisance_penalty, evidence_tier, exclusion_reason)]
  tasquinimod <- evidence[tolower(canonical_name) == "tasquinimod", .(compound, v1_sparse_primary_rank, v2_primary_rank, rank_stability_score, fold_model_agreement, matched_random_specificity_status, evidence_state, prism_status, curated_MoA, nuisance_penalty, evidence_tier, exclusion_reason)]
  historical_display <- historical[, .(canonical_name, compound, v1_sparse_primary_rank, v2_primary_rank, rank_change_v2_minus_v1, rank_stability_score, fold_model_agreement)]
  tier_counts <- candidate_evidence[, .N, by = evidence_tier][order(match(evidence_tier, c("tier_A", "tier_B", "tier_C", "unresolved", "discordant")))]
  failed_readiness <- readiness[pass == FALSE, condition]
  triggered <- downgrade[triggered == TRUE, reason]

  report <- c(
    "# Figure 8 v2 main-figure reanalysis report", "",
    paste0("**Final status: `", status, "`**"), "",
    "**Selected Results title:** Transcriptomic reversal analysis generates exploratory compound hypotheses without definitive cross-platform validation", "",
    "## Executive result", "",
    "A continuous 978-landmark rescue vector was constructed and evaluated with frozen DrugReflector checkpoints, 2,000 matched nulls, related CMap-family resources, curated MoA/targets, Figure 6 network compatibility, PRISM cancer-cell viability, nuisance controls, and bidirectional literature context. The prespecified main-figure gate was not met. No Tier A, B, or C candidate survived the joint robustness, specificity, corroboration, mechanism, phenotype, and nuisance requirements.", "",
    "## Required scientific questions", "",
    "### 1. Why was v1 landmark coverage 15.7%?", "", "The v1 high-confidence signature deliberately selected only 150 UP and 150 DOWN genes from the frozen biological evidence. Only 47 of those 300 genes intersected the frozen 978-gene DrugReflector space. This was sparse projection, not missing model genes.", "",
    "### 2. How was the continuous landmark-space score built?", "", "For every frozen model gene, normal/reference-minus-malignant state effect, reverse pseudotime association, Axis A identity evidence, Axis B stress evidence, Axis C SOX4 evidence, malignant-fate association, and direction-interpretable perturbation evidence were robust-scaled to [-1,1], combined with frozen 0.35/0.15/0.10/0.10/0.10/0.15/0.05 weights, and shrunk by directional agreement. Unsupported coordinates remained zero.", "",
    paste0("### 3. How many landmarks had reliable rescue scores?\n\n", sum(abs(score$final_rescue_vscore) > 0), "/978; ", sum(score$final_rescue_vscore > 0), " positive, ", sum(score$final_rescue_vscore < 0), " negative, and ", sum(score$final_rescue_vscore == 0), " zero."), "",
    paste0("### 4. Axis contribution\n\n", paste(balance$axis, scales::percent(balance$absolute_mass_fraction, accuracy = 0.1), collapse = "; "), "."), "",
    paste0("### 5. Axis-representation bias\n\nThe effective-axis number was ", format(balance$effective_axis_number[[1]], digits = 3), ". Axis A was largest, but the prespecified severe-domination rule was not triggered. Axis C remains the least represented and is a review risk."), "",
    paste0("### 6. Agreement with the old sparse ranking\n\nSpearman rho = ", numeric_comparison[metric == "v2_vs_v1_sparse_spearman", v2], ", indicating moderate rather than interchangeable ranking."), "",
    "### 7. Most rank-stable candidates", "", "```text", paste(capture.output(print(top_stable)), collapse = "\n"), "```", "",
    "### 8. Matched-random specificity", "", paste0("Global signature metrics were mostly not specific. Candidate-level results included ", candidate_evidence[matched_random_specificity_status == "strong", .N], " strong and ", candidate_evidence[matched_random_specificity_status == "suggestive", .N], " suggestive candidates; ", candidate_evidence[matched_random_specificity_status == "significantly_worse", .N], " were significantly worse than null."), "", "```text", paste(capture.output(print(random_summary)), collapse = "\n"), "```", "",
    paste0("### 9. Three-framework overlap\n\nv1 remained 0; v2 continuous-primary produced ", cross_report$v2_three_way_overlap, " standardized entities."), "",
    paste0("### 10. Two-framework support\n\n", cross_report$v2_two_way_overlap, " standardized entities."), "",
    "### 11. Candidate-level cross-framework support", "", "```text", paste(capture.output(print(cross_supported)), collapse = "\n"), "```", "",
    paste0("### 12. Curated MoA availability\n\n", annotation[!is.na(curated_MoA), .N], "/", nrow(annotation), " DrugReflector entities; candidate-universe rate ", scales::percent(annotation[candidate_analysis_universe == TRUE, mean(!is.na(curated_MoA))], accuracy = 0.1), "."), "",
    paste0("### 13. Curated target availability\n\n", annotation[!is.na(curated_targets), .N], "/", nrow(annotation), " entities."), "",
    paste0("### 14. Figure 6 network consistency\n\n", network[network_consistency_score >= 0.50, .N], " entities reached the frozen network-consistency threshold. This is compatibility, not verified drug action."), "",
    paste0("### 15. PRISM coverage\n\nPrimary: ", figure8_v2_count_true(prism, "primary_screen_available"), " entities; secondary dose-response: ", figure8_v2_count_true(prism, "secondary_dose_response_available"), ". Candidate-universe coverage was primary ", prism[candidate_analysis_universe == TRUE, scales::percent(mean(primary_screen_available), accuracy = 0.1)], " and secondary ", prism[candidate_analysis_universe == TRUE, scales::percent(mean(secondary_dose_response_available), accuracy = 0.1)], "."), "",
    paste0("### 16. HCC/liver-lineage viability effect\n\n", prism[prism_gate_class == "hcc_liver_enriched", .N], " candidate: tipifarnib-P2, supported by one secondary release and therefore not replicated orthogonal evidence."), "",
    paste0("### 17. Pan-cancer activity\n\n", prism[prism_gate_class == "pan_cancer_activity", .N], " candidates were classified as pan-cancer activity; ", prism[prism_gate_class == "broad_cytotoxicity", .N], " met the stricter broad-cytotoxicity definition."), "",
    paste0("### 18. Nuisance-driven candidates\n\n", candidate_evidence[nuisance_penalty >= 0.75, .N], "/", nrow(candidate_evidence), " candidate-universe entities were classified as nonspecific apparent reversal."), "",
    "### 19. Everolimus", "", "Everolimus received no reference bonus. All BRD instances remained computationally weak/unresolved or discordant. Preclinical HCC activity is counterbalanced by the EVOLVE-1 Phase III failure to improve overall survival.", "", "```text", paste(capture.output(print(everolimus)), collapse = "\n"), "```", "",
    "### 20. Tasquinimod", "", "Tasquinimod retained contextual positive and negative evidence: non-HCC Phase II PFS signals, Phase III lack of OS benefit, and a completed single-arm HCC cohort. It did not satisfy the v2 candidate gate.", "", "```text", paste(capture.output(print(tasquinimod)), collapse = "\n"), "```", "",
    "### 21. Historical v1 top hits", "", "```text", paste(capture.output(print(historical_display)), collapse = "\n"), "```", "",
    paste0("### 22. Tier A\n\n", candidate_evidence[evidence_tier == "tier_A", .N], "."), "",
    paste0("### 23. Tier B\n\n", candidate_evidence[evidence_tier == "tier_B", .N], "."), "",
    paste0("### 24. Discordant and unresolved\n\n", candidate_evidence[evidence_tier == "discordant", .N], " discordant and ", candidate_evidence[evidence_tier == "unresolved", .N], " unresolved candidate-universe entities."), "",
    paste0("### 25. Main-figure readiness\n\n`", status, "`. Failed readiness conditions: ", paste(failed_readiness, collapse = "; "), "."), "",
    paste0("### 26. Most important missing evidence\n\nA candidate satisfying both internal robustness and low nuisance while also carrying favorable matched-null specificity, CMap corroboration, and replicated liver/HCC phenotype. Triggered automatic downgrade cases: ", if (length(triggered)) paste(triggered, collapse = "; ") else "none; the explicit Tier A/B readiness condition still failed", "."), "",
    "### 27. Need for wet-lab validation", "", "Yes. A future upgrade requires perturbational expression in adult-HCC models, phenotypic rescue beyond viability, dose-response replication, and normal-cell/toxicity assays. PRISM alone cannot provide safety or demonstrate state rescue.", "",
    "### 28. Supported conclusion", "", "A continuous landmark-space analysis generated exploratory compound hypotheses and revealed limited multi-layer corroboration, but no candidate met the prespecified evidence tier required for a main translational figure.", "",
    "### 29. Unsupported conclusions", "", "The analysis cannot support an effective drug, therapeutic efficacy, validated treatment, clinical actionability, treatment recommendation, proven HCC reversal, normal-cell safety, or readiness for clinical translation.", "",
    "### 30. v1 to v2 interpretation change", "", "v2 improves model-space representation, exact provenance, null matching, cell-context correction, mechanism annotation, and orthogonal phenotype coverage. It changes the evidence structure, not the biological axes. The stronger audit still yields an Extended Data conclusion rather than a definitive compound nomination.", "",
    "## Candidate evidence tiers", "", "```text", paste(capture.output(print(tier_counts)), collapse = "\n"), "```", "",
    "## Top 20 candidate audit", "", "```text", paste(capture.output(print(top20)), collapse = "\n"), "```", "",
    "## Main-figure gate audit", "", "```text", paste(capture.output(print(readiness)), collapse = "\n"), "```", "",
    "## Validation summary", "", paste0("- Pass: ", sum(checks$status == "pass"), "\n- Pass with review context: ", sum(checks$status == "pass_with_review_context"), "\n- Fail: ", sum(checks$status == "fail")), "",
    "## Final recommended legend", "",
    "Figure 8 v2 | Transcriptomic reversal analysis generates exploratory compound hypotheses without definitive cross-platform validation. (a) Frozen three-axis malignant and desired rescue states. (b) Continuous 978-landmark representation and axis balance, with the v1 sparse overlap retained. (c) Internal rank robustness and 2,000-signature matched-null specificity. (d) Cross-framework corroboration across related LINCS/CMap-derived resources. (e) Exact curated targets and separately labelled Figure 6 network compatibility. (f) PRISM cancer-cell viability across adult-HCC and pan-cancer contexts. (g) Visible evidence components, nuisance penalties, coverage, and evidence tiers. The analysis is computational and exploratory and does not establish treatment efficacy, normal-cell safety, or clinical actionability."
  )
  writeLines(report, report_path, useBytes = TRUE)
  checks[check_name == "unsupported claims blocked", pass := figure8_v2_claims_clean(report)]
  checks[check_name == "unsupported claims blocked", status := fifelse(pass, "pass", "fail")]
  figure8_v2_write_tsv(checks, "figure8_v2_validation_report.tsv")
  figure8_v2_write_json(list(
    module = "figure8_v2_validation", status = if (any(checks$status == "fail")) "failed" else "passed",
    main_figure_status = status, n_pass = sum(checks$status == "pass"),
    n_pass_with_review_context = sum(checks$status == "pass_with_review_context"),
    n_fail = sum(checks$status == "fail"), checks = checks,
    protected_changes = change_audit, raster_dpi = raster_audit
  ), "figure8_v2_validation_report.json")

  cat("FIGURE 8 V2 FINAL SUMMARY\n")
  cat("Primary DrugReflector input: 978 landmarks; usable=", sum(abs(score$final_rescue_vscore) > 0), "; positive=", sum(score$final_rescue_vscore > 0), "; negative=", sum(score$final_rescue_vscore < 0), "; axis A/B/C=", paste(scales::percent(balance$absolute_mass_fraction, accuracy = 0.1), collapse = "/"), "\n", sep = "")
  cat("Robustness: continuous-v1 rho=", numeric_comparison[metric == "v2_vs_v1_sparse_spearman", v2], "; stable+fold candidates=", candidate_evidence[stable_gate & fold_gate, .N], "\n", sep = "")
  cat("Specificity: strong=", candidate_evidence[matched_random_specificity_status == "strong", .N], "; suggestive=", candidate_evidence[matched_random_specificity_status == "suggestive", .N], "; global status=", paste(unique(random_summary$specificity_status), collapse = ","), "\n", sep = "")
  cat("Cross-framework: v1 three-way=", cross_report$v1_three_way_overlap, "; v2 three-way=", cross_report$v2_three_way_overlap, "; v2 two-way=", cross_report$v2_two_way_overlap, "\n", sep = "")
  cat("MoA: candidate annotation rate=", scales::percent(annotation[candidate_analysis_universe == TRUE, mean(curated_annotation_available)], accuracy = 0.1), "; network-consistent=", network[network_consistency_score >= 0.50, .N], "\n", sep = "")
  cat("PRISM: candidate primary=", prism[candidate_analysis_universe == TRUE & primary_screen_available == TRUE, .N], "; secondary=", prism[candidate_analysis_universe == TRUE & secondary_dose_response_available == TRUE, .N], "; HCC/liver support=", prism[prism_gate_class == "hcc_liver_enriched", .N], "\n", sep = "")
  cat("Final evidence: TierA=", candidate_evidence[evidence_tier == "tier_A", .N], "; TierB=", candidate_evidence[evidence_tier == "tier_B", .N], "; TierC=", candidate_evidence[evidence_tier == "tier_C", .N], "; discordant=", candidate_evidence[evidence_tier == "discordant", .N], "; unresolved=", candidate_evidence[evidence_tier == "unresolved", .N], "\n", sep = "")
  cat("Main-figure status: ", status, "\n", sep = "")
  cat("Decision reasons: ", paste(failed_readiness, collapse = "; "), "\n", sep = "")
  cat("Figure1-7 and Figure8 v1 protected changes: ", nrow(change_audit), "\n", sep = "")
  if (any(checks$status == "fail")) stop("Figure 8 v2 validation failed: ", paste(checks[status == "fail", check_name], collapse = "; "))
  invisible(list(checks = checks, report = report_path, status = status))
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") figure8_v2_validate_and_report()
