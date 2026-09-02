root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_08_nuisance_literature.R")
if (!file.exists(script)) stop("Expected RED failure: nuisance/literature module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

stopifnot(figure8_v2_curated_nuisance_flag("protein synthesis inhibitor") == 1)
stopifnot(figure8_v2_curated_nuisance_flag("PPAR agonist") == 0)
stopifnot(figure8_v2_combined_nuisance_penalty(0.4, "protein synthesis inhibitor", "no_enriched_support") >= 0.9)
stopifnot(figure8_v2_combined_nuisance_penalty(0.4, "PPAR agonist", "broad_cytotoxicity") == 1)

stopifnot(figure8_v2_literature_direction("Phase III trial did not improve overall survival") == "negative")
stopifnot(figure8_v2_literature_direction("Compound inhibited HCC growth in xenografts") == "positive_preclinical")
stopifnot(figure8_v2_literature_direction("Pharmacokinetic study") == "contextual")
stopifnot(figure8_v2_literature_primary_score_weight() == 0)

candidate_set <- figure8_v2_literature_candidate_set(c("A", "B"), c("b", "C"), c("D"), c("A", "E"), c("F"))
stopifnot(identical(candidate_set, c("a", "b", "c", "d", "e", "f")))

cat("figure8_v2 nuisance/literature logic tests passed\n")
