suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/run_slingshot_module5_3.R <run_dir> [embedding=scanvi|hepatocyte_pca] [output_dir]", call. = FALSE)
}

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
embedding <- if (length(args) >= 2) args[[2]] else "scanvi"
output_dir <- if (length(args) >= 3) args[[3]] else file.path(run_dir, paste0("slingshot_", embedding))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(output_dir, "slingshot_status.json")
write_status <- function(status, message, extra = list()) {
  fields <- c(
    list(
      method = paste0("slingshot_", embedding),
      status = status,
      message = message,
      run_dir = run_dir
    ),
    extra
  )
  lines <- c("{")
  names_fields <- names(fields)
  for (i in seq_along(fields)) {
    value <- fields[[i]]
    if (is.numeric(value) || is.integer(value)) {
      rendered <- as.character(value)
    } else if (is.logical(value)) {
      rendered <- ifelse(value, "true", "false")
    } else {
      rendered <- paste0('"', gsub('"', '\\"', as.character(value)), '"')
    }
    comma <- if (i < length(fields)) "," else ""
    lines <- c(lines, paste0('  "', names_fields[[i]], '": ', rendered, comma))
  }
  lines <- c(lines, "}")
  writeLines(lines, status_path, useBytes = TRUE)
}

missing <- c()
for (pkg in c("SingleCellExperiment", "slingshot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}
if (length(missing) > 0) {
  write_status("missing_package", paste("Missing R package(s):", paste(missing, collapse = ", ")))
  stop(paste("Missing R package(s):", paste(missing, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

meta <- read_tsv_gz(file.path(run_dir, "cell_metadata.tsv.gz"))
embed_file <- switch(
  embedding,
  scanvi = "embedding_x_scanvi.tsv.gz",
  x_scanvi = "embedding_x_scanvi.tsv.gz",
  hepatocyte_pca = "embedding_hepatocyte_pca.tsv.gz",
  pca = "embedding_hepatocyte_pca.tsv.gz",
  stop("embedding must be one of: scanvi, x_scanvi, hepatocyte_pca, pca", call. = FALSE)
)
embed <- read_tsv_gz(file.path(run_dir, embed_file))
rownames(embed) <- embed$cell_id
embed_mat <- as.matrix(embed[meta$cell_id, setdiff(colnames(embed), "cell_id"), drop = FALSE])

cluster_labels <- as.character(meta$slingshot_cluster)
start_cluster <- unique(as.character(meta$slingshot_start_cluster))
start_cluster <- start_cluster[!is.na(start_cluster) & nzchar(start_cluster)][1]
end_clusters <- unique(unlist(strsplit(paste(unique(meta$slingshot_end_clusters), collapse = ","), ",")))
end_clusters <- end_clusters[!is.na(end_clusters) & nzchar(end_clusters)]

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Matrix::Matrix(0, nrow = 1, ncol = nrow(meta), sparse = TRUE)),
  colData = S4Vectors::DataFrame(meta)
)
colnames(sce) <- meta$cell_id
rownames(sce) <- "placeholder"
reducedDims(sce)[[embedding]] <- embed_mat

fit <- tryCatch(
  {
    slingshot::slingshot(
      sce,
      clusterLabels = cluster_labels,
      reducedDim = embedding,
      start.clus = start_cluster,
      end.clus = end_clusters,
      dist.method = "simple",
      shrink = FALSE,
      reweight = FALSE,
      reassign = FALSE,
      maxit = 10,
      approx_points = 150
    )
  },
  error = function(e) {
    write_status("failed", conditionMessage(e))
    stop(e)
  }
)
sce <- fit

pt <- as.data.frame(slingshot::slingPseudotime(sce))
pt$cell_id <- rownames(pt)
pt <- pt[, c("cell_id", setdiff(colnames(pt), "cell_id")), drop = FALSE]
write.table(pt, file.path(output_dir, "slingshot_pseudotime.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(sce, file.path(output_dir, "slingshot_sce.rds"))

write_status(
  "ok",
  "Slingshot completed.",
  list(
    n_cells = ncol(sce),
    start_cluster = start_cluster,
    end_clusters = paste(end_clusters, collapse = ","),
    n_lineages = ncol(pt) - 1
  )
)
