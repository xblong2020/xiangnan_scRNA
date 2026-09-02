root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_04_fetch_external_resources.R")
if (!file.exists(script)) stop("Expected RED failure: external-resource module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

plan <- figure8_v2_external_resource_plan()
stopifnot(nrow(plan) >= 11L)
stopifnot(all(grepl("^https://", plan$api_url)))
stopifnot(any(plan$filename == "Repurposing_Public_23Q2_LFC_COLLAPSED.csv"))
stopifnot(any(plan$filename == "secondary-screen-dose-response-curve-parameters.csv"))
stopifnot(any(plan$filename == "secondary-screen-replicate-collapsed-logfold-change.csv"))
stopifnot(any(plan$filename == "Model.csv" & plan$release == "DepMap 23Q2 Public v4"))
article_plan <- figure8_v2_plan_for_article(plan, 23600310L)
stopifnot(nrow(article_plan) == 6L, all(article_plan$article_id == 23600310L))
article <- figure8_v2_fetch_json("https://api.figshare.com/v2/articles/23600310")
stopifnot(article$id == 23600310L, article$version == 4L)

tmp <- tempfile(fileext = ".txt")
writeLines("figure8-v2", tmp, useBytes = TRUE)
expected <- unname(tools::md5sum(tmp))
record <- figure8_v2_validate_external_file(tmp, file.info(tmp)$size, expected)
stopifnot(record$status == "verified", record$observed_md5 == expected)
bad <- figure8_v2_validate_external_file(tmp, file.info(tmp)$size + 1, expected)
stopifnot(bad$status == "size_mismatch")
unlink(tmp)

empty_missing <- figure8_v2_missing_download_rows(data.table::data.table(
  resource_id = "x", release = "v1", filename = "x.tsv", role = "test", status = "verified"
))
stopifnot(nrow(empty_missing) == 0L)
stopifnot(all(c("resource_id", "release", "filename", "role", "status", "error") %in% names(empty_missing)))

cat("figure8_v2 external-resource logic tests passed\n")
