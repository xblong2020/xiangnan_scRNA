root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_02_analyze_drugreflector.R")
if (!file.exists(script)) stop("Expected RED failure: robustness module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

stopifnot(figure8_v2_rank_percentile(1, 9597) == 1)
stopifnot(figure8_v2_rank_percentile(9597, 9597) == 0)

stable <- figure8_v2_rank_stability(rep(10, 16), 9597, top_k = 200)
variable <- figure8_v2_rank_stability(c(1, 9597), 9597, top_k = 200)
stopifnot(isTRUE(all.equal(stable, 1)))
stopifnot(variable < 0.5)

stopifnot(isTRUE(all.equal(figure8_v2_fold_agreement(c(10, 10, 10), 9597), 1)))
stopifnot(figure8_v2_fold_agreement(c(1, 9597, 9597), 9597) == 0)

ordered <- figure8_v2_order_historical(
  data.table::data.table(canonical_name = c("cefepime", "dapivirine", "tasquinimod"), compound = c("C", "A", "T")),
  c("dapivirine", "cefepime", "tasquinimod")
)
stopifnot(identical(tolower(ordered$canonical_name), c("dapivirine", "cefepime", "tasquinimod")))

cat("figure8_v2 robustness logic tests passed\n")
