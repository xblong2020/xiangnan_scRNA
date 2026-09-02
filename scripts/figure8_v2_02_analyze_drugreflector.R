#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_rank_stability <- function(ranks, n_compounds, top_k = 200L) {
  ranks <- as.numeric(ranks)
  ranks <- ranks[is.finite(ranks)]
  if (!length(ranks)) return(NA_real_)
  iqr_component <- 1 - stats::IQR(ranks, type = 7) / max(1, n_compounds - 1)
  range_component <- 1 - (max(ranks) - min(ranks)) / max(1, n_compounds - 1)
  top_component <- mean(ranks <= top_k)
  pmax(0, pmin(1, 0.50 * iqr_component + 0.25 * range_component + 0.25 * top_component))
}

figure8_v2_fold_agreement <- function(fold_ranks, n_compounds) {
  fold_ranks <- as.numeric(fold_ranks)
  fold_ranks <- fold_ranks[is.finite(fold_ranks)]
  if (!length(fold_ranks)) return(NA_real_)
  pmax(0, pmin(1, 1 - (max(fold_ranks) - min(fold_ranks)) / max(1, n_compounds - 1)))
}

figure8_v2_order_historical <- function(x, name_order) {
  x <- as.data.table(copy(x))
  x[, historical_order__ := match(tolower(canonical_name), tolower(name_order))]
  setorder(x, historical_order__, compound, na.last = TRUE)
  x[, historical_order__ := NULL]
  x
}

figure8_v2_analyze_drugreflector <- function() {
  predictions <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_variant_predictions.tsv.gz"))
  folds <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_fold_predictions.tsv.gz"))
  manifest <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_signature_variant_manifest.tsv"))
  expected <- unique(manifest$signature_id)
  if (uniqueN(predictions$signature_id) != 16L || !setequal(unique(predictions$signature_id), expected)) stop("DrugReflector predictions do not contain all 16 frozen signatures")
  n_compounds <- uniqueN(predictions$compound)
  if (n_compounds != 9597L) stop("Unexpected DrugReflector compound universe")

  predictions[, rank_percentile := figure8_v2_rank_percentile(rank_1based, n_compounds)]
  primary_id <- "landmark_continuous_three_axis_rescue_vscore"
  primary <- predictions[signature_id == primary_id, .(
    compound, v2_primary_rank = rank_1based, v2_primary_rank_percentile = rank_percentile,
    v2_primary_probability = probability, v2_primary_logit = logit
  )]
  old_primary <- predictions[signature_id == "primary_three_axis", .(compound, v1_sparse_primary_rank = rank_1based, v1_sparse_primary_probability = probability)]

  stability <- predictions[, .(
    median_rank = median(rank_1based),
    rank_q1 = as.numeric(quantile(rank_1based, 0.25, type = 7)),
    rank_q3 = as.numeric(quantile(rank_1based, 0.75, type = 7)),
    rank_iqr = IQR(rank_1based, type = 7),
    best_rank = min(rank_1based),
    worst_rank = max(rank_1based),
    rank_range = max(rank_1based) - min(rank_1based),
    rank_sd = sd(rank_1based),
    top10_frequency = mean(rank_1based <= 10),
    top20_frequency = mean(rank_1based <= 20),
    top50_frequency = mean(rank_1based <= 50),
    top100_frequency = mean(rank_1based <= 100),
    top200_frequency = mean(rank_1based <= 200),
    rank_stability_score = figure8_v2_rank_stability(rank_1based, n_compounds, 200L),
    n_signature_versions = .N
  ), by = compound]

  primary_folds <- folds[signature_id == primary_id]
  fold_summary <- primary_folds[, .(
    fold_rank_min = min(fold_rank_1based),
    fold_rank_max = max(fold_rank_1based),
    fold_rank_iqr = IQR(fold_rank_1based, type = 7),
    fold_rank_sd = sd(fold_rank_1based),
    fold_model_agreement = figure8_v2_fold_agreement(fold_rank_1based, n_compounds)
  ), by = compound]
  full <- Reduce(function(x, y) merge(x, y, by = "compound", all = TRUE), list(stability, primary, old_primary, fold_summary))

  identity <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8_perturbagen_identity_map.tsv.gz"))
  identity <- identity[, .(compound, canonical_name, normalized_name, inchi_key, pubchem_cid, metadata_sources, metadata_conflict_flag)]
  full <- merge(full, identity, by = "compound", all.x = TRUE)
  full[is.na(canonical_name) | canonical_name == "", canonical_name := compound]
  full[, standardized_id := figure8_v2_entity_key(inchi_key, canonical_name, compound)]

  rank_wide <- dcast(predictions, compound ~ signature_id, value.var = "rank_1based")
  rank_names <- setdiff(names(rank_wide), "compound")
  setnames(rank_wide, rank_names, paste0("rank__", gsub("[^A-Za-z0-9]+", "_", rank_names)))
  full <- merge(full, rank_wide, by = "compound", all.x = TRUE)

  top200_union <- predictions[rank_1based <= 200, unique(compound)]
  historical_names <- c(
    "dapivirine", "tipiracil", "cefepime", "cisapride", "pomalidomide", "exalamide",
    "tg-100801", "tasquinimod", "levocetirizine", "escitalopram", "everolimus"
  )
  historical_ids <- full[tolower(canonical_name) %in% historical_names, unique(compound)]
  candidate_universe <- unique(c(top200_union, historical_ids))
  full[, `:=`(
    candidate_analysis_universe = compound %in% candidate_universe,
    historical_v1_top_or_reference = compound %in% historical_ids,
    candidate_origin = fifelse(compound %in% historical_ids, "historical_v1_or_reference", "top200_union_16_signatures")
  )]
  setorder(full, v2_primary_rank, compound)
  full[, v2_primary_priority_rank := .I]
  figure8_v2_write_tsv(full, "figure8_v2_drugreflector_full_ranking.tsv.gz", compress = TRUE)
  figure8_v2_write_tsv(full[candidate_analysis_universe == TRUE, .(compound)], "figure8_v2_candidate_watchlist.tsv")

  correlations <- merge(
    predictions[signature_id == primary_id, .(compound, primary_rank = rank_1based)],
    predictions[signature_id != primary_id, .(compound, signature_id, comparison_rank = rank_1based)],
    by = "compound", allow.cartesian = TRUE
  )[, .(
    n_compounds = .N,
    spearman_rho = cor(primary_rank, comparison_rank, method = "spearman"),
    kendall_tau = cor(primary_rank, comparison_rank, method = "kendall"),
    top100_overlap = length(intersect(compound[primary_rank <= 100], compound[comparison_rank <= 100])),
    top200_overlap = length(intersect(compound[primary_rank <= 200], compound[comparison_rank <= 200]))
  ), by = signature_id]
  figure8_v2_write_tsv(correlations, "figure8_v2_signature_rank_correlations.tsv")

  historical <- full[tolower(canonical_name) %in% historical_names, .(
    canonical_name, compound, inchi_key, v1_sparse_primary_rank, v2_primary_rank,
    rank_change_v2_minus_v1 = v2_primary_rank - v1_sparse_primary_rank,
    median_rank, rank_q1, rank_q3, best_rank, worst_rank, rank_stability_score, fold_model_agreement
  )]
  historical <- figure8_v2_order_historical(historical, historical_names)
  figure8_v2_write_tsv(historical, "figure8_v2_historical_top_hits_rank_audit.tsv")

  primary_old_rho <- correlations[signature_id == "primary_three_axis", spearman_rho]
  figure8_v2_write_json(list(
    module = "figure8_v2_drugreflector_robustness", status = "completed",
    seed = FIGURE8_V2_SEED, n_compounds = n_compounds, n_signatures = uniqueN(predictions$signature_id),
    primary_signature = primary_id, old_primary_spearman_rho = primary_old_rho,
    candidate_universe_n = length(candidate_universe),
    stability_formula = "0.50*(1-IQR/(N-1)) + 0.25*(1-range/(N-1)) + 0.25*top200_frequency",
    fold_agreement_formula = "1 - fold rank range/(N-1)",
    interpretation_boundary = "Internal model robustness only; no efficacy or phenotype claim."
  ), "figure8_v2_drugreflector_robustness_report.json")
  invisible(full)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_analyze_drugreflector()
  cat("FIGURE8_V2_DRUGREFLECTOR compounds=", nrow(result), " candidates=", sum(result$candidate_analysis_universe), "\n", sep = "")
}
