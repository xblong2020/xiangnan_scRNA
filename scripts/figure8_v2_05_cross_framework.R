#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_cross_state <- function(dr, l1000, clue, l1000_profiled, clue_profiled, opposition) {
  if (isTRUE(opposition)) return("discordant")
  support_n <- sum(c(isTRUE(dr), isTRUE(l1000), isTRUE(clue)))
  if (support_n == 3L) return("three_framework_support")
  if (support_n == 2L) return("two_framework_support")
  if (isTRUE(dr) && !isTRUE(l1000_profiled) && !isTRUE(clue_profiled)) return("drugreflector_only_external_not_available")
  if (isTRUE(dr)) return("drugreflector_only")
  if (isTRUE(l1000) || isTRUE(clue)) return("external_only")
  if (isTRUE(l1000_profiled) || isTRUE(clue_profiled)) return("profiled_no_directional_support")
  "not_available"
}

figure8_v2_cross_framework <- function() {
  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  v1 <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8e_method_support_matrix.tsv"))

  dr <- ranking[, .(
    canonical_name = canonical_name[which.max(!is.na(canonical_name) & canonical_name != "")][1],
    inchi_key = inchi_key[which.max(!is.na(inchi_key) & inchi_key != "")][1],
    brd_ids = paste(sort(unique(compound)), collapse = ";"),
    v2_drugreflector_rank = min(v2_primary_rank),
    v2_drugreflector_support = min(v2_primary_rank) <= 200,
    v2_drugreflector_profiled = TRUE,
    identity_mapping_conflict = any(metadata_conflict_flag, na.rm = TRUE)
  ), by = standardized_id]
  external <- v1[, .(
    standardized_id, external_canonical_name = canonical_name, external_inchi_key = inchi_key,
    mapping_level, mapping_confidence, external_mapping_conflict = mapping_conflict,
    L1000FWD = as.logical(L1000FWD), CLUE = as.logical(CLUE),
    l1000_result_available = as.logical(l1000_result_listed), clue_result_available = as.logical(clue_result_available),
    l1000_rank, clue_tau, external_strong_opposition = as.logical(any_strong_opposition),
    v1_drugreflector_support = as.logical(DrugReflector), v1_drugreflector_rank = drugreflector_rank,
    v1_liver_context_score = liver_context_score, v1_liver_support = as.logical(liver_support)
  )]
  cross <- merge(external, dr, by = "standardized_id", all = TRUE)
  for (nm in c("L1000FWD", "CLUE", "l1000_result_available", "clue_result_available", "external_strong_opposition", "v1_drugreflector_support", "v2_drugreflector_support", "v2_drugreflector_profiled", "identity_mapping_conflict", "external_mapping_conflict")) {
    cross[is.na(get(nm)), (nm) := FALSE]
  }
  cross[, canonical_name := fifelse(!is.na(canonical_name) & canonical_name != "", canonical_name, external_canonical_name)]
  cross[, inchi_key := fifelse(!is.na(inchi_key) & inchi_key != "", inchi_key, external_inchi_key)]
  cross[, mapping_conflict := identity_mapping_conflict | external_mapping_conflict]
  cross[, evidence_state := mapply(
    figure8_v2_cross_state,
    v2_drugreflector_support, L1000FWD, CLUE,
    l1000_result_available, clue_result_available, external_strong_opposition,
    USE.NAMES = FALSE
  )]
  cross[, n_support_frameworks := as.integer(v2_drugreflector_support) + as.integer(L1000FWD) + as.integer(CLUE)]
  cross[, `:=`(
    l1000_status = fifelse(!l1000_result_available, "not_available", fifelse(L1000FWD, "same_direction_support", fifelse(external_strong_opposition, "opposition", "profiled_no_support"))),
    clue_status = fifelse(!clue_result_available, "not_available", fifelse(CLUE, "same_direction_support", fifelse(external_strong_opposition, "opposition", "profiled_no_support"))),
    approximate_name_matching_used = FALSE,
    framework_relationship = "cross-framework corroboration across LINCS/CMap-derived resources"
  )]
  setorder(cross, -n_support_frameworks, v2_drugreflector_rank, standardized_id)
  figure8_v2_write_tsv(cross, "figure8_v2_cross_framework_concordance.tsv")
  conflicts <- cross[mapping_conflict == TRUE | is.na(standardized_id), .(
    standardized_id, canonical_name, inchi_key, brd_ids, mapping_level, mapping_confidence,
    identity_mapping_conflict, external_mapping_conflict, conflict_reason = fifelse(is.na(standardized_id), "no_standardized_identifier", "source_identifier_conflict")
  )]
  figure8_v2_write_tsv(conflicts, "figure8_v2_compound_mapping_conflicts.tsv")

  levels <- rbindlist(list(
    data.table(mapping_level_audit = "standardized_entity", v2_three_way = cross[v2_drugreflector_support & L1000FWD & CLUE, .N], v2_two_way = cross[n_support_frameworks == 2, .N]),
    data.table(mapping_level_audit = "exact_BRD_instance", v2_three_way = NA_integer_, v2_two_way = NA_integer_),
    data.table(mapping_level_audit = "canonical_name_or_InChIKey", v2_three_way = cross[v2_drugreflector_support & L1000FWD & CLUE & grepl("^(INCHI|NAME):", standardized_id), .N], v2_two_way = cross[n_support_frameworks == 2 & grepl("^(INCHI|NAME):", standardized_id), .N])
  ))
  figure8_v2_write_tsv(levels, "figure8_v2_cross_framework_mapping_level_summary.tsv")

  v1_three_way <- v1[DrugReflector & L1000FWD & CLUE, .N]
  v2_three_way <- cross[v2_drugreflector_support & L1000FWD & CLUE, .N]
  v2_two_way <- cross[n_support_frameworks == 2, .N]
  figure8_v2_write_json(list(
    module = "figure8_v2_cross_framework", status = "completed",
    direction_rules = list(
      DrugReflector = "v2 continuous-primary rank <=200 is discovery support",
      L1000FWD = "frozen similar-to-rescue direction after v1 query-orientation audit",
      CLUE = "frozen positive tau support after v1 desired-rescue query-orientation audit"
    ),
    v1_three_way_overlap = v1_three_way,
    v2_three_way_overlap = v2_three_way,
    v2_two_way_overlap = v2_two_way,
    mapping_conflicts = nrow(conflicts),
    approximate_name_matching_used = FALSE,
    interpretation_boundary = "L1000FWD and CLUE are related LINCS/CMap-derived corroboration resources, not independent validation platforms."
  ), "figure8_v2_cross_framework_report.json")
  invisible(cross)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_cross_framework()
  cat("FIGURE8_V2_CROSS entities=", nrow(result), " two_way=", sum(result$n_support_frameworks == 2), " three_way=", sum(result$n_support_frameworks == 3), "\n", sep = "")
}

