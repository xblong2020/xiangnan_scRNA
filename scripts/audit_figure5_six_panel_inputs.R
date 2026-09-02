#!/usr/bin/env Rscript

# Capture a read-only checksum baseline for the existing corrected Figure 5
# artefacts before the six-panel namespace is populated.  This manifest is
# intentionally limited to data/metadata/figure outputs; source code may gain
# new helpers while the protected analytical results remain byte-identical.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."),
                winslash = "/", mustWork = TRUE)
} else normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(root, "scripts", "figure5_temporal_core.R"))
old <- figure5_onset_fix_paths(root, create = FALSE)
new <- figure5_six_panel_paths(root, create = TRUE)

protected_roots <- c(old$metadata, old$processed, old$figures, old$preview)
files <- unique(unlist(lapply(protected_roots, function(path) {
  if (!dir.exists(path)) return(character())
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
})))
files <- files[file.info(files)$isdir %in% FALSE]
if (!length(files)) stop("No corrected Figure 5 artefacts found to protect", call. = FALSE)

manifest <- rbindlist(lapply(files, function(path) {
  info <- file.info(path)
  data.table(
    file_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = unname(tools::md5sum(path))
  )
}), fill = TRUE)
setorder(manifest, file_path)
out <- file.path(new$metadata, "figure5_six_panel_protected_original_manifest.tsv")
fwrite(manifest, out, sep = "\t", quote = FALSE, na = "")
cat(sprintf("Protected %d corrected artefacts in %s\n", nrow(manifest), out))

