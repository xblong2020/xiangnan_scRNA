root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_09_integrate_evidence.R")
if (!file.exists(script)) stop("Expected RED failure: integration module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

stopifnot(figure8_v2_rank_percentile(1, 9597) == 1)
stopifnot(figure8_v2_rank_percentile(9597, 9597) == 0)

stopifnot(figure8_v2_specificity_score("strong") == 1)
stopifnot(figure8_v2_specificity_score("suggestive") == 0.6)
stopifnot(figure8_v2_specificity_score("not_specific") == 0)
stopifnot(figure8_v2_specificity_score("significantly_worse") == 0)

components <- c(DR_score = 1, robustness_score = 0.8, signature_specificity_score = 1, cross_framework_score = NA, network_moa_score = 0.5, prism_phenotype_score = NA)
summary <- figure8_v2_integrated_summary(components, nuisance_penalty = 0.2)
stopifnot(summary$coverage_confidence == 0.75)
stopifnot(is.finite(summary$conservative_score), is.finite(summary$coverage_aware_score))

row <- data.table::as.data.table(as.list(components))
extracted <- figure8_v2_component_vector(row, names(components))
stopifnot(identical(names(extracted), names(components)))
stopifnot(isTRUE(all.equal(unname(extracted), unname(components))))

stopifnot(figure8_v2_tier_order("tier_A") < figure8_v2_tier_order("tier_B"))
stopifnot(figure8_v2_tier_order("tier_B") < figure8_v2_tier_order("tier_C"))
stopifnot(figure8_v2_tier_order("unresolved") < figure8_v2_tier_order("discordant"))
stopifnot(figure8_v2_reference_bonus("everolimus") == 0)

cat("figure8_v2 integration logic tests passed\n")
