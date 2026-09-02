root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_03_matched_random.R")
if (!file.exists(script)) stop("Expected RED failure: matched-random module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

toy <- data.table::data.table(
  gene = paste0("G", 1:8),
  final_rescue_vscore = c(0.8, 0.4, -0.9, -0.3, 0.2, -0.5, 0, 0),
  primary_axis_stratum = c("A", "A", "B", "B", "C", "C", "zero", "zero"),
  mean_expression = 1:8,
  expression_variance = seq(0.1, 0.8, 0.1),
  detection_rate = seq(0.2, 0.9, 0.1)
)
set.seed(123)
null <- figure8_v2_permute_matched(toy)
stopifnot(identical(null$gene, toy$gene))
stopifnot(sum(null$vscore > 0) == sum(toy$final_rescue_vscore > 0))
stopifnot(sum(null$vscore < 0) == sum(toy$final_rescue_vscore < 0))
stopifnot(sum(null$vscore == 0) == sum(toy$final_rescue_vscore == 0))
for (stratum in c("A", "B", "C")) {
  observed <- sort(abs(toy[primary_axis_stratum == stratum, final_rescue_vscore]))
  randomized <- sort(abs(null[primary_axis_stratum == stratum, vscore]))
  stopifnot(identical(observed, randomized))
}

qc <- figure8_v2_matching_qc(toy, null)
stopifnot(qc$pass)
stopifnot(qc$nonzero_count_difference == 0)
stopifnot(qc$positive_count_difference == 0)
stopifnot(qc$negative_count_difference == 0)
stopifnot(qc$absolute_weight_ks <= 0.10)
stopifnot(max(abs(c(qc$expression_smd, qc$variance_smd, qc$detection_smd))) <= 0.10)

p <- figure8_v2_empirical_p(0.9, c(0.1, 0.2, 0.3, 0.4))
stopifnot(p$p_two_sided == 0.4)
stopifnot(figure8_v2_directional_specificity_label(0.01, 0.99, 0.02) == "strong")
stopifnot(figure8_v2_directional_specificity_label(0.99, 0.01, 0.02) == "significantly_worse")
stopifnot(figure8_v2_directional_specificity_label(0.30, 0.70, 0.60) == "not_specific")

top_metrics <- figure8_v2_random_top_metrics(
  data.table::data.table(signature_id = c("N1", "N1", "N2"), compound = c("A", "B", "C"), rank_1based = c(1L, 2L, 1L)),
  old_top200 = c("A"), historical = c("B")
)
stopifnot(is.double(top_metrics$best_historical_rank))
stopifnot(top_metrics[signature_id == "N1", best_historical_rank] == 2)
stopifnot(is.na(top_metrics[signature_id == "N2", best_historical_rank]))

cat("figure8_v2 matched-random logic tests passed\n")
