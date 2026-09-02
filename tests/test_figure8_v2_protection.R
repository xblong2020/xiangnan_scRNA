`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1]])) y else x

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(project_root, "AGENTS.md")))

script_path <- file.path(project_root, "scripts", "figure8_v2_00_freeze_audit.R")
if (!file.exists(script_path)) stop("Expected RED failure: figure8_v2_00_freeze_audit.R is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script_path, local = FALSE)

stopifnot(exists("figure8_v2_protected_files"))
stopifnot(exists("figure8_v2_hash_manifest"))
stopifnot(exists("figure8_v2_design_changes"))

paths <- figure8_v2_protected_files(project_root)
rel <- gsub("\\\\", "/", paths)

stopifnot(length(paths) > 0L)
stopifnot(any(grepl("^scripts/figure[1-7]", rel, ignore.case = TRUE)))
stopifnot(any(grepl("^metadata/driver/figure8_transcriptomic_reversal/", rel)))
stopifnot(any(grepl("trajectory|celloracle|scenic|sctenifold", rel, ignore.case = TRUE)))
stopifnot(!any(grepl("figure8_transcriptomic_reversal_v2_mainfigure|figure8_v2_", rel, ignore.case = TRUE)))

sample_paths <- paths[seq_len(min(3L, length(paths)))]
manifest <- figure8_v2_hash_manifest(sample_paths, project_root, workers = 1L)
stopifnot(identical(names(manifest), c("file_path", "size_bytes", "modified_utc", "md5")))
stopifnot(nrow(manifest) == length(sample_paths))
stopifnot(all(nchar(manifest$md5) == 32L))
stopifnot(all(file.exists(file.path(project_root, manifest$file_path))))

changes <- figure8_v2_design_changes()
stopifnot(all(c("component", "v1", "v2", "reason") %in% names(changes)))
stopifnot(any(changes$component == "matched-random specificity" & grepl("1000", changes$v1)))
stopifnot(any(changes$component == "DrugReflector input" & grepl("continuous", changes$v2)))

cat("figure8_v2 protection tests passed\n")
