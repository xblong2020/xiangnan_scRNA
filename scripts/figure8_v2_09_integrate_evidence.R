#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_positive_weights <- function() c(
  DR_score = 0.20,
  robustness_score = 0.20,
  signature_specificity_score = 0.20,
  cross_framework_score = 0.15,
  network_moa_score = 0.15,
  prism_phenotype_score = 0.10
)

figure8_v2_specificity_score <- function(status) {
  status <- as.character(status)
  ifelse(status == "strong", 1, ifelse(status == "suggestive", 0.6, ifelse(status %in% c("not_specific", "significantly_worse"), 0, NA_real_)))
}

figure8_v2_integrated_summary <- function(components, nuisance_penalty) {
  weights <- figure8_v2_positive_weights()
  components <- components[names(weights)]
  evidence <- figure8_v2_evidence_summary(components, weights)
  penalty <- if (is.finite(nuisance_penalty)) 0.25 * nuisance_penalty else 0
  list(
    conservative_score = pmax(0, evidence$conservative_score - penalty),
    coverage_aware_score = pmax(0, evidence$coverage_aware_score - penalty),
    coverage_confidence = evidence$coverage_confidence,
    missing_dimensions = paste(names(weights)[evidence$missing_dimensions], collapse = ";")
  )
}

figure8_v2_tier_order <- function(tier) {
  match(as.character(tier), c("tier_A", "tier_B", "tier_C", "unresolved", "discordant"))
}

figure8_v2_reference_bonus <- function(candidate_name) 0

figure8_v2_component_vector <- function(row, component_names = names(figure8_v2_positive_weights())) {
  values <- as.numeric(unlist(as.list(row[, ..component_names]), use.names = FALSE))
  stats::setNames(values, component_names)
}

figure8_v2_collapse_literature <- function(literature) {
  literature <- as.data.table(literature)
  literature[, candidate_key := figure8_v2_safe_name(candidate)]
  collapse_titles <- function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    if (length(x)) paste(x, collapse = " | ") else NA_character_
  }
  literature[, .(
    positive_literature = collapse_titles(title[grepl("^positive", evidence_direction)]),
    negative_literature = collapse_titles(title[grepl("negative|inconclusive", evidence_direction)]),
    clinical_evidence = collapse_titles(title[grepl("phase|clinical|trial", tolower(evidence_type)) | !is.na(NCT)]),
    literature_rows = .N,
    positive_literature_count = sum(grepl("^positive", evidence_direction)),
    negative_literature_count = sum(grepl("negative|inconclusive", evidence_direction))
  ), by = candidate_key]
}

figure8_v2_integrate_evidence <- function() {
  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  specificity <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_candidate_matched_random_specificity.tsv"))
  cross <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_concordance.tsv"))
  annotation <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_compound_moa_target_annotation.tsv"))
  network <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_network_consistency.tsv"))
  prism <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_prism_viability.tsv"))
  nuisance <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_nuisance_penalties.tsv"))
  literature <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_candidate_literature_evidence.tsv"))

  x <- merge(ranking, specificity, by.x = "compound", by.y = "compound", all.x = TRUE)
  x <- merge(x, cross[, .(
    standardized_id, L1000FWD, CLUE, l1000_result_available, clue_result_available,
    l1000_status, clue_status, n_support_frameworks, evidence_state,
    external_strong_opposition, v1_liver_support, v1_liver_context_score
  )], by = "standardized_id", all.x = TRUE)
  x <- merge(x, annotation[, .(
    BRD_ID, curated_MoA, curated_targets, clinical_status, approved_status,
    annotation_confidence, mapping_conflict, curated_annotation_available
  )], by.x = "compound", by.y = "BRD_ID", all.x = TRUE)
  x <- merge(x, network[, .(
    BRD_ID, network_consistency_score, network_evidence_coverage,
    direct_network_targets, one_step_targets, pathway_targets,
    compatible_axes, network_consistent_inferred_mechanism
  )], by.x = "compound", by.y = "BRD_ID", all.x = TRUE)
  x <- merge(x, prism[, .(
    BRD_ID, prism_status, primary_screen_available, secondary_dose_response_available,
    prism_phenotype_class, prism_gate_class, prism_phenotype_score,
    primary_adult_hcc_median_sensitivity, primary_pan_cancer_median_sensitivity,
    primary_liver_minus_non_liver, n_primary_adult_hcc_lines,
    secondary_adult_hcc_median_auc, secondary_pan_cancer_median_auc,
    secondary_liver_minus_non_liver, n_secondary_adult_hcc_lines,
    normal_cell_safety_established
  )], by.x = "compound", by.y = "BRD_ID", all.x = TRUE)
  x <- merge(x, nuisance[, .(
    BRD_ID, proliferation_penalty, generic_stress_penalty, dna_damage_penalty,
    translation_inhibition_penalty, mitochondrial_toxicity_penalty,
    pan_cytotoxicity_penalty, model_nuisance_penalty, curated_nuisance_flag,
    nuisance_penalty, nuisance_interpretation
  )], by.x = "compound", by.y = "BRD_ID", all.x = TRUE)

  x[, `:=`(
    DR_score = figure8_v2_rank_percentile(v2_primary_rank, 9597),
    robustness_score = rowMeans(cbind(rank_stability_score, fold_model_agreement), na.rm = TRUE),
    signature_specificity_score = figure8_v2_specificity_score(matched_random_specificity_status)
  )]
  x[!is.finite(robustness_score), robustness_score := NA_real_]
  x[, external_profile_count := as.integer(l1000_result_available %in% TRUE) + as.integer(clue_result_available %in% TRUE)]
  x[, external_support_count := as.integer(L1000FWD %in% TRUE) + as.integer(CLUE %in% TRUE)]
  x[, cross_framework_score := fifelse(external_profile_count > 0, external_support_count / external_profile_count, NA_real_)]
  x[, network_moa_score := fifelse(
    curated_annotation_available %in% TRUE,
    pmin(1, 0.5 + 0.5 * fifelse(is.finite(network_consistency_score), network_consistency_score, 0)),
    NA_real_
  )]
  x[, prism_phenotype_score := as.numeric(prism_phenotype_score)]

  component_cols <- names(figure8_v2_positive_weights())
  summaries <- lapply(seq_len(nrow(x)), function(idx) {
    figure8_v2_integrated_summary(figure8_v2_component_vector(x[idx], component_cols), x$nuisance_penalty[idx])
  })
  x[, `:=`(
    conservative_score = vapply(summaries, `[[`, numeric(1), "conservative_score"),
    coverage_aware_score = vapply(summaries, `[[`, numeric(1), "coverage_aware_score"),
    evidence_coverage = vapply(summaries, `[[`, numeric(1), "coverage_confidence"),
    missing_positive_dimensions = vapply(summaries, `[[`, character(1), "missing_dimensions")
  )]

  x[, stable_gate := is.finite(rank_stability_score) & rank_stability_score >= 0.75]
  x[, fold_gate := is.finite(fold_model_agreement) & fold_model_agreement >= 0.75]
  x[, specificity_gate := matched_random_specificity_status %in% c("strong", "suggestive")]
  x[, cmap_gate := external_support_count >= 1]
  x[, mechanism_gate := is.finite(network_moa_score) & network_moa_score >= 0.50]
  x[, nuisance_gate := !is.finite(nuisance_penalty) | nuisance_penalty < 0.75]
  x[, prism_gate := is.na(prism_gate_class) | !prism_gate_class %in% c("broad_cytotoxicity", "discordant_phenotype")]
  x[, n_passed_core_gates := rowSums(cbind(stable_gate, fold_gate, specificity_gate, cmap_gate, mechanism_gate, nuisance_gate, prism_gate), na.rm = TRUE)]
  tier_rows <- lapply(seq_len(nrow(x)), function(idx) {
    figure8_v2_assign_tier(
      rank_stability = x$rank_stability_score[idx], fold_agreement = x$fold_model_agreement[idx],
      specificity_p = if (x$specificity_gate[idx]) x$matched_random_specificity_p[idx] else NA_real_,
      cmap_support = x$cmap_gate[idx], cmap_profiled = x$external_profile_count[idx] > 0,
      network_moa_score = x$network_moa_score[idx], nuisance_penalty = x$nuisance_penalty[idx],
      prism_class = x$prism_gate_class[idx], coverage = x$evidence_coverage[idx],
      strong_opposition = x$external_strong_opposition[idx]
    )
  })
  x[, `:=`(
    evidence_tier = vapply(tier_rows, `[[`, character(1), "tier"),
    exclusion_reason = vapply(tier_rows, `[[`, character(1), "failed_gates"),
    reference_bonus = vapply(canonical_name, figure8_v2_reference_bonus, numeric(1))
  )]
  if (any(x$reference_bonus != 0)) stop("Reference compounds received a forbidden score bonus")
  x[, tier_order := figure8_v2_tier_order(evidence_tier)]
  setorder(x, tier_order, -n_passed_core_gates, -conservative_score, -coverage_aware_score, v2_primary_rank, compound)
  x[, candidate_priority_rank := .I]
  figure8_v2_write_tsv(x, "figure8_v2_integrated_candidate_evidence.tsv")

  lit_summary <- figure8_v2_collapse_literature(literature)
  top20 <- x[candidate_analysis_universe == TRUE][1:20]
  top20[, candidate_key := figure8_v2_safe_name(canonical_name)]
  top20 <- merge(top20, lit_summary, by = "candidate_key", all.x = TRUE)
  top20_audit <- top20[, .(
    candidate_priority_rank, BRD_ID = compound, canonical_name, InChIKey = inchi_key,
    DrugReflector_median_rank = median_rank, rank_IQR = rank_iqr,
    rank_range = rank_range, v2_primary_rank, model_agreement = fold_model_agreement,
    matched_random_specificity = matched_random_specificity_status,
    matched_random_empirical_P = matched_random_specificity_p,
    L1000FWD_status = l1000_status, CLUE_status = clue_status,
    PRISM_status = prism_status, liver_HCC_phenotype = prism_phenotype_class,
    pan_cancer_phenotype = fifelse(prism_gate_class == "broad_cytotoxicity", "broad_cytotoxicity", fifelse(prism_gate_class == "pan_cancer_activity", "pan_cancer_activity", "not_broad_or_unavailable")),
    MoA = curated_MoA, targets = curated_targets,
    Figure6_network_compatibility = network_consistent_inferred_mechanism,
    network_consistency_score, nuisance_penalty,
    positive_literature, negative_literature, clinical_evidence,
    evidence_tier, exclusion_reason, evidence_coverage,
    normal_cell_safety_established = FALSE
  )]
  figure8_v2_write_tsv(top20_audit, "figure8_v2_top20_candidate_audit.tsv")

  v1 <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8h_candidate_ranking_full.tsv"))
  v1_top <- v1[order(final_priority_rank)]
  v2_top <- x[candidate_analysis_universe == TRUE]
  rank_corr <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_signature_rank_correlations.tsv"))
  axis_balance <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_landmark_axis_balance.tsv"))
  random_summary <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_random_specificity_summary.tsv"))
  numeric_comparison <- data.table(
    metric = c(
      "primary_input_dimension", "usable_landmark_representation", "top10_candidate_overlap", "top50_candidate_overlap",
      "top100_candidate_overlap", "top200_candidate_overlap", "v2_vs_v1_sparse_spearman",
      "median_rank_stability", "median_model_agreement", "matched_random_global_specificity",
      "cross_framework_three_way_count", "cross_framework_two_way_count", "PRISM_HCC_liver_support_count",
      "MoA_annotation_rate_candidate_universe", "discordant_count", "unresolved_count", "Tier_A_count", "Tier_B_count", "Tier_C_count"
    ),
    v1 = c(
      "300 sparse genes", "47/300 (15.7%)", NA, NA, NA, "127 primary/sensitivity overlap", NA,
      median(v1$top200_frequency, na.rm = TRUE), median(v1$model_agreement, na.rm = TRUE), "specificity gate FALSE",
      0, 14, 0, 0, 1041, 8226, 0, 0, 0
    ),
    v2 = c(
      "978 model landmarks", "839/978 non-zero continuous scores", length(intersect(v1_top$compound[1:10], v2_top$compound[1:10])), length(intersect(v1_top$compound[1:50], v2_top$compound[1:50])),
      length(intersect(v1_top$compound[1:100], v2_top$compound[1:100])), length(intersect(v1_top$compound[1:200], v2_top$compound[1:200])), rank_corr[signature_id == "primary_three_axis", spearman_rho],
      median(x$rank_stability_score, na.rm = TRUE), median(x$fold_model_agreement, na.rm = TRUE), paste(unique(random_summary$specificity_status), collapse = ";"),
      sum(x$n_support_frameworks == 3, na.rm = TRUE), sum(x$n_support_frameworks == 2, na.rm = TRUE), sum(x$prism_gate_class == "hcc_liver_enriched", na.rm = TRUE),
      mean(x[candidate_analysis_universe == TRUE, curated_annotation_available], na.rm = TRUE), sum(x$evidence_tier == "discordant"), sum(x$evidence_tier == "unresolved"), sum(x$evidence_tier == "tier_A"), sum(x$evidence_tier == "tier_B"), sum(x$evidence_tier == "tier_C")
    )
  )
  figure8_v2_write_tsv(numeric_comparison, "figure8_v1_vs_v2_numeric_comparison.tsv")

  changes <- c(
    "# Figure 8 v1 to v2 interpretation changes", "",
    "- v1 sparse landmark coverage remains a historical limitation; v2 uses the actual 978-coordinate model space with unsupported coordinates fixed at zero.",
    paste0("- Continuous-vs-v1 sparse rank Spearman rho: ", format(rank_corr[signature_id == "primary_three_axis", spearman_rho], digits = 4), "."),
    paste0("- CMap-family three-way overlap changed from 0 to ", sum(x$n_support_frameworks == 3, na.rm = TRUE), "; these frameworks remain related and are not called independent validation."),
    paste0("- PRISM found ", sum(x$prism_gate_class == "hcc_liver_enriched", na.rm = TRUE), " HCC/liver-enriched candidate(s); normal-cell safety remains unknown."),
    paste0("- Candidate-level matched-null specificity: ", sum(x$matched_random_specificity_status == "strong", na.rm = TRUE), " strong and ", sum(x$matched_random_specificity_status == "suggestive", na.rm = TRUE), " suggestive; global signature metrics remain mostly non-specific."),
    "- Literature is contextual only and contributes zero primary score weight."
  )
  writeLines(changes, file.path(FIGURE8_V2_METADATA, "figure8_v1_vs_v2_interpretation_changes.md"), useBytes = TRUE)

  candidate_x <- x[candidate_analysis_universe == TRUE]
  tier_counts <- candidate_x[, .N, by = evidence_tier]
  readiness <- data.table(
    condition_id = 1:10,
    condition = c(
      "continuous landmark-space input established", "no severe single-axis domination", "internal rank robustness interpretable",
      "matched-random specificity has partial favorable support", "three-way overlap is not the sole validation gate",
      "candidate-level CMap-family corroboration exists", "major-candidate MoA/targets partially available",
      "orthogonal PRISM evidence exists", "top candidates are not all nuisance dominated", "at least one Tier A or multiple Tier B"
    ),
    pass = c(
      file.exists(file.path(FIGURE8_V2_METADATA, "figure8_v2_gene_level_rescue_vscore.tsv")),
      !axis_balance$severe_axis_domination[[1]],
      sum(candidate_x$stable_gate & candidate_x$fold_gate, na.rm = TRUE) > 0,
      sum(candidate_x$specificity_gate, na.rm = TRUE) > 0,
      TRUE,
      sum(candidate_x$cmap_gate, na.rm = TRUE) > 0,
      mean(candidate_x$curated_annotation_available, na.rm = TRUE) >= 0.30,
      sum(candidate_x$primary_screen_available | candidate_x$secondary_dose_response_available, na.rm = TRUE) > 0,
      mean(top20$nuisance_gate, na.rm = TRUE) > 0,
      sum(candidate_x$evidence_tier == "tier_A") >= 1 | sum(candidate_x$evidence_tier == "tier_B") >= 2
    )
  )
  automatic_downgrade <- data.table(
    case = 1:6,
    reason = c("matched-random specificity completely failed", "continuous and old rankings almost completely inconsistent", "no orthogonal phenotype support", "all top hits generic nuisance/cytotoxicity", "MoA/targets mostly unavailable", "only DrugReflector supports candidates"),
    triggered = c(
      sum(candidate_x$specificity_gate, na.rm = TRUE) == 0,
      rank_corr[signature_id == "primary_three_axis", spearman_rho] < 0.20,
      sum(candidate_x$primary_screen_available | candidate_x$secondary_dose_response_available, na.rm = TRUE) == 0,
      all(!top20$nuisance_gate),
      mean(candidate_x$curated_annotation_available, na.rm = TRUE) < 0.50,
      sum(candidate_x$cmap_gate, na.rm = TRUE) == 0
    )
  )
  main_status <- if (all(readiness$pass) && !any(automatic_downgrade$triggered)) "MAIN_FIGURE_READY" else "EXTENDED_DATA_ONLY"
  readiness[, main_figure_status := main_status]
  figure8_v2_write_tsv(readiness, "figure8_v2_main_figure_readiness_gates.tsv")
  figure8_v2_write_tsv(automatic_downgrade, "figure8_v2_automatic_downgrade_audit.tsv")
  figure8_v2_write_json(list(
    module = "figure8_v2_integrated_evidence", status = "completed",
    main_figure_status = main_status, tier_counts = tier_counts,
    readiness = readiness, automatic_downgrade = automatic_downgrade,
    score_definition = "reports/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_integrated_score_definition.md",
    reference_bonus = 0, literature_primary_score_weight = 0,
    interpretation_boundary = "Exploratory computational and cancer-cell phenotype evidence; no effective treatment or clinical recommendation claim."
  ), "figure8_v2_integrated_evidence_report.json")
  invisible(list(evidence = x, top20 = top20_audit, readiness = readiness, status = main_status))
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_integrate_evidence()
  cat("FIGURE8_V2_INTEGRATED status=", result$status, " TierA=", sum(result$evidence$evidence_tier == "tier_A"), " TierB=", sum(result$evidence$evidence_tier == "tier_B"), " TierC=", sum(result$evidence$evidence_tier == "tier_C"), "\n", sep = "")
}
