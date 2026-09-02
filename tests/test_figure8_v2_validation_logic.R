root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_12_validate_report.R")
if (!file.exists(script)) stop("Expected RED failure: validation/report module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

checks <- figure8_v2_required_check_names()
stopifnot(length(checks) >= 30L)
stopifnot(all(c("Figure1-7 protected", "Figure8 v1 protected", "978 landmarks loaded", "random benchmark >=1000", "unsupported claims blocked") %in% checks))

before <- data.table::data.table(file_path = c("a", "b"), size_bytes = c(1, 2), md5 = c("x", "y"))
after <- data.table::data.table(file_path = c("a", "b"), size_bytes = c(1, 3), md5 = c("x", "z"))
audit <- figure8_v2_compare_hashes(before, after)
stopifnot(audit[file_path == "a", status] == "unchanged")
stopifnot(audit[file_path == "b", status] == "changed")

stopifnot(!figure8_v2_claims_clean("This is an effective drug"))
stopifnot(figure8_v2_claims_clean("This is an exploratory predicted candidate hypothesis"))
stopifnot(figure8_v2_claims_clean("The analysis cannot support a validated treatment or a clinically actionable compound."))
stopifnot(figure8_v2_count_true(data.table::data.table(flag = c(TRUE, FALSE, NA)), "flag") == 1L)

main_png <- file.path("figures", "driver", "figure8_transcriptomic_reversal_v2_mainfigure", "figure8_v2_mainfigure_a_to_g.png")
if (file.exists(main_png)) {
  info <- figure8_v2_image_info_ascii(main_png)
  stopifnot(info$width[[1]] > 0, info$height[[1]] > 0)
}

cat("figure8_v2 validation logic tests passed\n")
