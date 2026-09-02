root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_06_moa_network.R")
if (!file.exists(script)) stop("Expected RED failure: MoA/network module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

ids <- figure8_v2_brd_core(c("BRD-K71847383-001-12-5", "BRD-A00077618-236-07-6", NA))
stopifnot(identical(ids, c("BRD-K71847383", "BRD-A00077618", NA_character_)))

score <- figure8_v2_network_score(direct_target = 1, one_step = 1, pathway_overlap = 0.5, axis_compatibility = 1)
stopifnot(isTRUE(all.equal(score, 0.90)))

compat <- figure8_v2_axis_compatibility("PPARA", "PPAR agonist")
stopifnot(compat$score == 1, compat$axis == "axis_A_identity")
compat_bad <- figure8_v2_axis_compatibility("HNF4A", "HNF4A inhibitor")
stopifnot(compat_bad$score == 0)

edges <- figure8_v2_target_edges(data.table::data.table(
  BRD_ID = c("B1", "B2"), canonical_name = c("drug1", "drug2"),
  curated_targets = c("HNF4A, PPARA", NA), source = c("official", NA),
  annotation_confidence = c("high", NA)
))
stopifnot(nrow(edges) == 2L, setequal(edges$target, c("HNF4A", "PPARA")))

cat("figure8_v2 MoA/network logic tests passed\n")
