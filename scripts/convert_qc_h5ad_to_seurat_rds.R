root <- normalizePath(getwd(), mustWork = TRUE)
qc_h5ad_root <- file.path(root, "data", "processed", "qc_h5ad")
out_root <- file.path(root, "data", "processed", "qc_seurat_rds")
meta_root <- file.path(root, "metadata", "qc")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

files <- list.files(qc_h5ad_root, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  dataset = character(),
  source_h5ad = character(),
  output_rds = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

for (source in files) {
  dataset <- basename(dirname(source))
  label <- sub("\\.h5ad$", "", basename(source))
  out_dir <- file.path(out_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  h5s <- file.path(out_dir, paste0(label, ".h5Seurat"))
  rds <- file.path(out_dir, paste0(label, ".rds"))
  if (file.exists(rds)) {
    manifest <- rbind(manifest, data.frame(dataset, source, rds, "skipped_existing", "", stringsAsFactors = FALSE))
    message("SKIP ", rds)
    next
  }
  status <- "complete"
  msg <- ""
  tryCatch({
    message("CONVERT ", source)
    Convert(source, dest = "h5seurat", overwrite = TRUE, filename = h5s, verbose = FALSE)
    obj <- LoadH5Seurat(h5s, verbose = FALSE)
    saveRDS(obj, rds, compress = "gzip")
    message("WROTE ", rds)
  }, error = function(e) {
    status <<- "failed"
    msg <<- conditionMessage(e)
    message("FAILED ", source, ": ", msg)
  })
  manifest <- rbind(manifest, data.frame(dataset, source, rds, status, msg, stringsAsFactors = FALSE))
}

write.table(
  manifest,
  file.path(meta_root, "qc_seurat_rds_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("WROTE manifest")
