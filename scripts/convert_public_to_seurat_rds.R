root <- normalizePath(getwd(), mustWork = TRUE)
out_root <- file.path(root, "data", "processed", "seurat_rds")
meta_root <- file.path(root, "metadata", "conversion")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

manifest <- data.frame(
  dataset = character(),
  source = character(),
  output = character(),
  class = character(),
  cells = numeric(),
  features = numeric(),
  status = character(),
  stringsAsFactors = FALSE
)

add_manifest <- function(dataset, source, output, object, status = "complete") {
  dims <- tryCatch(dim(object), error = function(e) c(NA, NA))
  manifest <<- rbind(
    manifest,
    data.frame(
      dataset = dataset,
      source = source,
      output = output,
      class = paste(class(object), collapse = "|"),
      cells = ifelse(length(dims) >= 2, dims[[2]], NA),
      features = ifelse(length(dims) >= 1, dims[[1]], NA),
      status = status,
      stringsAsFactors = FALSE
    )
  )
}

save_existing_rds <- function(dataset, source, output_name) {
  out_dir <- file.path(out_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, output_name)
  if (file.exists(out)) {
    obj <- readRDS(out)
    add_manifest(dataset, source, out, obj, "skipped_existing")
    message("SKIP ", out)
    return(invisible(out))
  }
  message("READ ", source)
  if (grepl("\\.gz$", source, ignore.case = TRUE)) {
    con <- gzcon(gzfile(source, "rb"))
    on.exit(close(con), add = TRUE)
    obj <- readRDS(con)
  } else {
    obj <- readRDS(source)
  }
  saveRDS(obj, out, compress = "gzip")
  add_manifest(dataset, source, out, obj)
  message("WROTE ", out)
  invisible(out)
}

read_10x_mtx_as_seurat <- function(dataset, matrix_path, features_path, barcodes_path, output_name) {
  out_dir <- file.path(out_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, output_name)
  if (file.exists(out)) {
    obj <- readRDS(out)
    add_manifest(dataset, matrix_path, out, obj, "skipped_existing")
    message("SKIP ", out)
    return(invisible(out))
  }
  message("READ 10x ", dataset)
  mat <- readMM(matrix_path)
  features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  gene_names <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  rownames(mat) <- make.unique(as.character(gene_names))
  colnames(mat) <- as.character(barcodes[[1]])
  cell_sums <- Matrix::colSums(mat)
  mat <- mat[, cell_sums > 0, drop = FALSE]
  colnames(mat) <- paste(dataset, colnames(mat), sep = "_")
  obj <- CreateSeuratObject(counts = mat, project = dataset)
  obj$dataset <- dataset
  saveRDS(obj, out, compress = "gzip")
  add_manifest(dataset, matrix_path, out, obj)
  message("WROTE ", out)
  invisible(out)
}

public <- file.path(root, "data", "public")

save_existing_rds(
  "GSE202379",
  file.path(public, "geo", "GSE202379", "GSE202379_SeuratObject_AllCells.rds.gz"),
  "GSE202379_SeuratObject_AllCells.rds"
)

save_existing_rds(
  "GSE174748",
  file.path(public, "geo", "GSE174748", "GSE174748_hl_nuclei.rds.gz"),
  "GSE174748_hl_nuclei.rds"
)

hcc_dir <- file.path(public, "figshare", "HCC_atlas")
for (source in list.files(hcc_dir, pattern = "\\.rds$", full.names = TRUE)) {
  save_existing_rds("HCC_atlas", source, basename(source))
}

read_10x_mtx_as_seurat(
  "GSE151530",
  file.path(public, "geo", "GSE151530", "GSE151530_matrix.mtx.gz"),
  file.path(public, "geo", "GSE151530", "GSE151530_genes.tsv.gz"),
  file.path(public, "geo", "GSE151530", "GSE151530_barcodes.tsv.gz"),
  "GSE151530_raw_counts.rds"
)

write.table(
  manifest,
  file.path(meta_root, "seurat_rds_conversion_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("WROTE manifest")
