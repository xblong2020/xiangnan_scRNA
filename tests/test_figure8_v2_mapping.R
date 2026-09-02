root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_05_cross_framework.R")
if (!file.exists(script)) stop("Expected RED failure: cross-framework module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

stopifnot(figure8_v2_cross_state(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) == "three_framework_support")
stopifnot(figure8_v2_cross_state(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) == "drugreflector_only_external_not_available")
stopifnot(figure8_v2_cross_state(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) == "profiled_no_directional_support")
stopifnot(figure8_v2_cross_state(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE) == "discordant")

script_text <- paste(readLines(script, warn = FALSE), collapse = "\n")
stopifnot(!grepl("adist|agrep|stringdist|fuzzy", script_text, ignore.case = TRUE))

cat("figure8_v2 mapping logic tests passed\n")
