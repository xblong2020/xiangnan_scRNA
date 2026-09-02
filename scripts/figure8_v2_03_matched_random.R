#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_primary_axis_stratum <- function(x) {
  x <- as.data.table(x)
  axis_matrix <- abs(as.matrix(x[, .(axis_A_component, axis_B_component, axis_C_component)]))
  axis_names <- c("axis_A", "axis_B", "axis_C")
  maximum <- apply(axis_matrix, 1, max)
  winner <- axis_names[max.col(axis_matrix, ties.method = "first")]
  winner[maximum <= 1e-12 & abs(x$malignant_state_component) > 1e-12] <- "malignant_state"
  winner[maximum <= 1e-12 & abs(x$malignant_state_component) <= 1e-12 & abs(x$final_rescue_vscore) > 0] <- "state_trajectory_only"
  winner[abs(x$final_rescue_vscore) <= 1e-12] <- "zero"
  paste(winner, fifelse(x$final_rescue_vscore > 0, "positive", fifelse(x$final_rescue_vscore < 0, "negative", "zero")), sep = "__")
}

figure8_v2_permute_matched <- function(observed) {
  observed <- as.data.table(copy(observed))
  if (!"primary_axis_stratum" %in% names(observed)) observed[, primary_axis_stratum := figure8_v2_primary_axis_stratum(observed)]
  observed[, vscore := sample(final_rescue_vscore, .N, replace = FALSE), by = primary_axis_stratum]
  observed[, .(gene, primary_axis_stratum, vscore, mean_expression, expression_variance, detection_rate)]
}

figure8_v2_smd <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  pooled <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (!is.finite(pooled) || pooled <= .Machine$double.eps) return(ifelse(isTRUE(all.equal(mean(a), mean(b))), 0, Inf))
  (mean(b) - mean(a)) / pooled
}

figure8_v2_matching_qc <- function(observed, randomized) {
  observed <- as.data.table(observed)
  randomized <- as.data.table(randomized)
  obs_nonzero <- observed[abs(final_rescue_vscore) > 1e-12]
  rnd_nonzero <- randomized[abs(vscore) > 1e-12]
  ks <- suppressWarnings(stats::ks.test(abs(obs_nonzero$final_rescue_vscore), abs(rnd_nonzero$vscore))$statistic[[1]])
  expression_smd <- figure8_v2_smd(obs_nonzero$mean_expression, rnd_nonzero$mean_expression)
  variance_smd <- figure8_v2_smd(obs_nonzero$expression_variance, rnd_nonzero$expression_variance)
  detection_smd <- figure8_v2_smd(obs_nonzero$detection_rate, rnd_nonzero$detection_rate)
  axis_observed <- table(observed$primary_axis_stratum)
  axis_random <- table(randomized$primary_axis_stratum)
  axis_names <- union(names(axis_observed), names(axis_random))
  axis_diff <- max(abs(axis_observed[axis_names] - axis_random[axis_names]), na.rm = TRUE)
  row <- data.table(
    nonzero_count_difference = nrow(rnd_nonzero) - nrow(obs_nonzero),
    positive_count_difference = sum(randomized$vscore > 0) - sum(observed$final_rescue_vscore > 0),
    negative_count_difference = sum(randomized$vscore < 0) - sum(observed$final_rescue_vscore < 0),
    axis_stratum_max_count_difference = axis_diff,
    absolute_weight_ks = ks,
    expression_smd = expression_smd,
    variance_smd = variance_smd,
    detection_smd = detection_smd
  )
  row[, pass := nonzero_count_difference == 0 & positive_count_difference == 0 & negative_count_difference == 0 &
    axis_stratum_max_count_difference == 0 & absolute_weight_ks <= 0.10 &
    max(abs(c(expression_smd, variance_smd, detection_smd))) <= 0.10]
  row
}

figure8_v2_random_top_metrics <- function(top, old_top200, historical) {
  as.data.table(top)[, .(
    top200_overlap_v1 = as.integer(sum(compound %in% old_top200)),
    best_historical_rank = {
      values <- as.numeric(rank_1based[compound %in% historical])
      if (length(values)) min(values) else NA_real_
    }
  ), by = signature_id]
}

figure8_v2_prepare_matched_random <- function(n_random = 2000L, seed = FIGURE8_V2_SEED) {
  set.seed(seed)
  score <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_gene_level_rescue_vscore.tsv"))
  model <- score[order(model_gene_order), .(
    gene, final_rescue_vscore, axis_A_component, axis_B_component, axis_C_component,
    malignant_state_component, mean_expression, expression_variance, detection_rate
  )]
  if (nrow(model) != 978L) stop("Matched-null source must contain exactly 978 landmarks")
  model[, primary_axis_stratum := figure8_v2_primary_axis_stratum(model)]
  matrix <- matrix(0, nrow = n_random, ncol = nrow(model), dimnames = list(sprintf("matched_null_%04d", seq_len(n_random)), model$gene))
  qc_rows <- vector("list", n_random)
  for (idx in seq_len(n_random)) {
    randomized <- figure8_v2_permute_matched(model)
    matrix[idx, ] <- randomized$vscore
    qc <- figure8_v2_matching_qc(model, randomized)
    qc[, signature_id := rownames(matrix)[idx]]
    qc_rows[[idx]] <- qc
  }
  qc <- rbindlist(qc_rows)
  if (!all(qc$pass)) stop("At least one matched-null signature violates frozen tolerances")
  hashes <- apply(matrix, 1, function(x) digest::digest(x, algo = "md5", serialize = TRUE))
  if (uniqueN(hashes) < n_random * 0.99) stop("Matched-null permutations are insufficiently unique")
  wide <- as.data.table(matrix, keep.rownames = "signature_id")
  figure8_v2_write_tsv(wide, "figure8_v2_matched_random_signatures_wide.tsv.gz", directory = FIGURE8_V2_DATA, compress = TRUE)
  qc[, `:=`(
    random_seed = seed,
    vector_md5 = hashes,
    nonzero_landmarks = sum(abs(model$final_rescue_vscore) > 0),
    positive_landmarks = sum(model$final_rescue_vscore > 0),
    negative_landmarks = sum(model$final_rescue_vscore < 0)
  )]
  figure8_v2_write_tsv(qc, "figure8_v2_matched_random_manifest.tsv.gz", compress = TRUE)
  invisible(qc)
}

figure8_v2_summarize_matched_random <- function() {
  inference <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_inference_summary.tsv.gz"))
  top <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_top_predictions.tsv.gz"))
  watched <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_watchlist_predictions.tsv.gz"))
  primary_pred <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_variant_predictions.tsv.gz"))[signature_id == "landmark_continuous_three_axis_rescue_vscore"]
  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  old_top200 <- ranking[v1_sparse_primary_rank <= 200, compound]
  historical <- ranking[historical_v1_top_or_reference == TRUE, compound]

  random_metrics <- merge(inference, figure8_v2_random_top_metrics(top, old_top200, historical), by = "signature_id", all.x = TRUE)

  observed_order <- primary_pred[order(rank_1based)]
  top_compound <- observed_order$compound[[1]]
  observed <- list(
    max_probability = max(primary_pred$probability),
    top10_probability_concentration = sum(observed_order[1:10, probability]),
    top100_probability_concentration = sum(observed_order[1:100, probability]),
    top200_probability_concentration = sum(observed_order[1:200, probability]),
    top200_overlap_v1 = sum(observed_order[1:200, compound] %in% old_top200),
    best_historical_rank = min(primary_pred[compound %in% historical, rank_1based]),
    top_candidate_model_agreement = ranking[compound == top_compound, fold_model_agreement]
  )
  metrics <- names(observed)
  summary <- rbindlist(lapply(metrics, function(metric) {
    null <- random_metrics[[metric]]
    value <- observed[[metric]]
    transformed_value <- if (metric == "best_historical_rank") -value else value
    transformed_null <- if (metric == "best_historical_rank") -null else null
    p <- figure8_v2_empirical_p(transformed_value, transformed_null)
    data.table(
      metric = metric, observed_value = value, random_median = median(null, na.rm = TRUE),
      random_q05 = as.numeric(quantile(null, 0.05, na.rm = TRUE)), random_q95 = as.numeric(quantile(null, 0.95, na.rm = TRUE)),
      empirical_p_upper = p$p_upper, empirical_p_lower = p$p_lower,
      empirical_p_two_sided = p$p_two_sided,
      specificity_status = figure8_v2_directional_specificity_label(p$p_upper, p$p_lower, p$p_two_sided)
    )
  }))
  figure8_v2_write_tsv(summary, "figure8_v2_random_specificity_summary.tsv")

  observed_candidate <- primary_pred[, .(
    compound, observed_rank = rank_1based,
    observed_rank_percentile = 1 - (rank_1based - 1) / (uniqueN(primary_pred$compound) - 1),
    observed_probability = probability
  )]
  watched[, null_rank_percentile := 1 - (rank_1based - 1) / (uniqueN(primary_pred$compound) - 1)]
  candidate_specificity <- watched[, {
    observed_row <- observed_candidate[compound == .BY$compound]
    p_rank <- figure8_v2_empirical_p(observed_row$observed_rank_percentile, null_rank_percentile)
    p_prob <- figure8_v2_empirical_p(observed_row$observed_probability, probability)
    .(
      observed_rank = observed_row$observed_rank,
      observed_rank_percentile = observed_row$observed_rank_percentile,
      observed_probability = observed_row$observed_probability,
      null_median_rank = median(rank_1based),
      null_rank_p_upper = p_rank$p_upper,
      null_rank_p_lower = p_rank$p_lower,
      null_rank_p_two_sided = p_rank$p_two_sided,
      null_probability_p_upper = p_prob$p_upper,
      null_probability_p_lower = p_prob$p_lower,
      null_probability_p_two_sided = p_prob$p_two_sided,
      matched_random_specificity_p = max(p_rank$p_two_sided, p_prob$p_two_sided),
      matched_random_specificity_status = fifelse(
        p_rank$p_upper <= 0.10 & p_prob$p_upper <= 0.10 & max(p_rank$p_two_sided, p_prob$p_two_sided) <= 0.10,
        fifelse(max(p_rank$p_two_sided, p_prob$p_two_sided) < 0.05, "strong", "suggestive"),
        fifelse(p_rank$p_lower < 0.05 | p_prob$p_lower < 0.05, "significantly_worse", "not_specific")
      )
    )
  }, by = compound]
  figure8_v2_write_tsv(candidate_specificity, "figure8_v2_candidate_matched_random_specificity.tsv")

  random_metrics[, `:=`(
    cmap_family_support_count = NA_integer_,
    prism_support_count = NA_integer_,
    response_class_concentration = NA_real_,
    integrated_evidence_score = NA_real_
  )]
  figure8_v2_write_tsv(random_metrics, "figure8_v2_matched_random_benchmark.tsv.gz", compress = TRUE)
  figure8_v2_write_json(list(
    module = "figure8_v2_matched_random", status = "completed_pending_external_metric_enrichment",
    n_random = nrow(inference), seed = FIGURE8_V2_SEED,
    matching = "same 978-gene universe; exact nonzero/sign/axis-stratum counts; within-stratum vscore permutation",
    empirical_p = "two-sided add-one formula",
    summary = summary
  ), "figure8_v2_matched_random_report.json")
  invisible(summary)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  stage_arg <- grep("^--stage=", args, value = TRUE)
  stage <- if (length(stage_arg)) sub("^--stage=", "", stage_arg[[1]]) else "prepare"
  if (stage == "prepare") {
    result <- figure8_v2_prepare_matched_random()
    cat("FIGURE8_V2_RANDOM_PREPARED n=", nrow(result), " pass=", sum(result$pass), "\n", sep = "")
  } else if (stage == "summarize") {
    result <- figure8_v2_summarize_matched_random()
    print(result)
  } else stop("Unknown stage: ", stage)
}
