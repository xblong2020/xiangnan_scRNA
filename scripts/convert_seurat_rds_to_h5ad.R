root <- normalizePath(getwd(), mustWork = TRUE)
seurat_root <- file.path(root, "data", "processed", "seurat_rds")
h5ad_root <- file.path(root, "data", "processed", "h5ad_from_seurat")
meta_root <- file.path(root, "metadata", "conversion")
dir.create(h5ad_root, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

rds_files <- list.files(seurat_root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  source = character(),
  h5seurat = character(),
  h5ad = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

add_manifest <- function(source, h5seurat, h5ad, status, message = "") {
  manifest <<- rbind(
    manifest,
    data.frame(
      source = source,
      h5seurat = h5seurat,
      h5ad = h5ad,
      status = status,
      message = message,
      stringsAsFactors = FALSE
    )
  )
}

for (source in rds_files) {
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(seurat_root)), "/?"), "", gsub("\\\\", "/", normalizePath(source)))
  rel_no_ext <- sub("\\.rds$", "", rel)
  h5seurat <- file.path(h5ad_root, paste0(rel_no_ext, ".h5Seurat"))
  h5ad <- file.path(h5ad_root, paste0(rel_no_ext, ".h5ad"))
  dir.create(dirname(h5ad), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(h5ad)) {
    add_manifest(source, h5seurat, h5ad, "skipped_existing", "")
    message("SKIP ", h5ad)
    next
  }

  status <- "complete"
  msg <- ""
  tryCatch({
    message("READ ", source)
    obj <- readRDS(source)
    if (grepl("GSE202379", source)) {
      obj@reductions <- list()
      obj@graphs <- list()
      DefaultAssay(obj) <- "RNA"
      if ("SCT" %in% Assays(obj)) {
        obj[["SCT"]] <- NULL
      }
    }
    obj <- UpdateSeuratObject(obj)
    DefaultAssay(obj) <- DefaultAssay(obj)
    SaveH5Seurat(obj, filename = h5seurat, overwrite = TRUE, verbose = FALSE)
    Convert(h5seurat, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
    message("WROTE ", h5ad)
  }, error = function(e) {
    status <<- "failed"
    msg <<- conditionMessage(e)
    message("FAILED ", source, ": ", msg)
  })
  add_manifest(source, h5seurat, h5ad, status, msg)
}

write.table(
  manifest,
  file.path(meta_root, "seurat_to_h5ad_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("WROTE manifest")
