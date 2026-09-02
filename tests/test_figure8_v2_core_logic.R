root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
common <- file.path(root, "scripts", "figure8_v2_common.R")
if (!file.exists(common)) stop("Expected RED failure: figure8_v2_common.R is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(common, local = FALSE)

v1_library <- figure8_v2_existing_r_library()
stopifnot(dir.exists(v1_library))

scaled <- figure8_v2_robust_unit(c(-100, -1, 0, 1, 100))
stopifnot(length(scaled) == 5L, all(is.finite(scaled)))
stopifnot(min(scaled) >= -1, max(scaled) <= 1, scaled[[3]] == 0)

agree <- figure8_v2_agreement_shrunk(c(1, -1), c(0.5, 0.5))
stopifnot(isTRUE(all.equal(agree$weighted_mean, 0)))
stopifnot(isTRUE(all.equal(agree$agreement, 0)))
stopifnot(isTRUE(all.equal(agree$score, 0)))

concordant <- figure8_v2_agreement_shrunk(c(1, 0.5), c(0.5, 0.5))
stopifnot(concordant$score > 0.5, concordant$agreement == 1)

keys <- figure8_v2_entity_key(
  inchi_key = c("AAAA-BBBB", NA, NA),
  canonical_name = c("ignored", "Drug A", NA),
  brd_id = c("BRD-1", "BRD-2", "BRD-3")
)
stopifnot(identical(keys, c("INCHI:AAAA-BBBB", "NAME:druga", "BRD:BRD-3")))
placeholder_key <- figure8_v2_entity_key("restricted", "QL-X-138", "BRD-U33728988")
stopifnot(placeholder_key == "NAME:qlx138")

p <- figure8_v2_empirical_p(observed = 5, null = 1:4)
stopifnot(isTRUE(all.equal(p$p_upper, 0.2)))
stopifnot(isTRUE(all.equal(p$p_lower, 1.0)))
stopifnot(isTRUE(all.equal(p$p_two_sided, 0.4)))
stopifnot(figure8_v2_specificity_label(0.049) == "strong")
stopifnot(figure8_v2_specificity_label(0.08) == "suggestive")
stopifnot(figure8_v2_specificity_label(0.11) == "not_specific")

evidence <- figure8_v2_evidence_summary(c(1, NA, 0.5), c(0.4, 0.3, 0.3))
stopifnot(isTRUE(all.equal(evidence$conservative_score, 0.55)))
stopifnot(isTRUE(all.equal(evidence$coverage_aware_score, 0.55 / 0.7)))
stopifnot(isTRUE(all.equal(evidence$coverage_confidence, 0.7)))

stopifnot(figure8_v2_prism_class(0.90, 0.50, 0.20, 5) == "hcc_liver_enriched")
stopifnot(figure8_v2_prism_class(0.90, 0.90, 0.01, 5) == "broad_cytotoxicity")
stopifnot(figure8_v2_prism_class(0.90, 0.50, 0.20, 2) == "insufficient_liver_lines")
stopifnot(is.na(figure8_v2_prism_score(NA_real_, NA_real_, NA_real_)))

tier_a <- figure8_v2_assign_tier(
  rank_stability = 0.80, fold_agreement = 0.80, specificity_p = 0.04,
  cmap_support = TRUE, cmap_profiled = TRUE, network_moa_score = 0.60,
  nuisance_penalty = 0.20, prism_class = "hcc_liver_enriched", coverage = 0.85,
  strong_opposition = FALSE
)
stopifnot(tier_a$tier == "tier_A")

tier_b <- figure8_v2_assign_tier(
  rank_stability = 0.80, fold_agreement = 0.80, specificity_p = 0.08,
  cmap_support = TRUE, cmap_profiled = TRUE, network_moa_score = 0.60,
  nuisance_penalty = 0.20, prism_class = NA_character_, coverage = 0.70,
  strong_opposition = FALSE
)
stopifnot(tier_b$tier == "tier_B")

discordant <- figure8_v2_assign_tier(
  rank_stability = 0.90, fold_agreement = 0.90, specificity_p = 0.01,
  cmap_support = TRUE, cmap_profiled = TRUE, network_moa_score = 0.90,
  nuisance_penalty = 0.90, prism_class = "broad_cytotoxicity", coverage = 1,
  strong_opposition = FALSE
)
stopifnot(discordant$tier == "discordant")

cat("figure8_v2 core logic tests passed\n")
