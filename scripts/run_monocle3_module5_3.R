suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/run_monocle3_module5_3.R <run_dir> [output_dir]", call. = FALSE)
}

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- if (length(args) >= 2) args[[2]] else file.path(run_dir, "monocle3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(output_dir, "monocle3_status.json")
write_status <- function(status, message, extra = list()) {
  fields <- c(
    list(
      method = "monocle3_learn_graph_order_cells",
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

if (!requireNamespace("monocle3", quietly = TRUE)) {
  write_status("missing_package", "R package monocle3 is not installed.")
  stop("R package monocle3 is not installed.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(monocle3)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

meta <- read_tsv_gz(file.path(run_dir, "cell_metadata.tsv.gz"))
genes <- read.delim(file.path(run_dir, "genes.tsv"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
counts <- readMM(gzfile(file.path(run_dir, "counts_gene_by_cell.mtx.gz")))
rownames(counts) <- genes$gene_id
colnames(counts) <- meta$cell_id

gene_meta <- data.frame(gene_short_name = genes$gene_short_name, row.names = genes$gene_id)
cell_meta <- meta
rownames(cell_meta) <- cell_meta$cell_id
cell_meta <- cell_meta[colnames(counts), , drop = FALSE]

cds <- monocle3::new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

umap <- read_tsv_gz(file.path(run_dir, "embedding_umap.tsv.gz"))
rownames(umap) <- umap$cell_id
umap <- as.matrix(umap[colnames(cds), setdiff(colnames(umap), "cell_id"), drop = FALSE])
colnames(umap) <- c("UMAP_1", "UMAP_2")
reducedDims(cds)$UMAP <- umap

root_cells <- cell_meta$cell_id[cell_meta$monocle3_root_cell %in% c(TRUE, "TRUE", "True", "true", 1, "1")]
if (length(root_cells) == 0) {
  write_status("failed", "No monocle3 root cells found in cell_metadata.tsv.gz.")
  stop("No monocle3 root cells found.", call. = FALSE)
}

cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", k = 30)
cds <- monocle3::learn_graph(cds, use_partition = FALSE)
cds <- monocle3::order_cells(cds, root_cells = root_cells)

pt <- monocle3::pseudotime(cds)
out <- data.frame(
  cell_id = names(pt),
  monocle3_pseudotime = as.numeric(pt),
  monocle3_is_finite = is.finite(as.numeric(pt)),
  stringsAsFactors = FALSE
)
cell_meta_out <- cell_meta[, c("cell_id", "trajectory_root_end_role", "cell_disease_stage", "sample_disease_stage", "slingshot_cluster"), drop = FALSE]
out <- merge(
  out,
  cell_meta_out,
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
write.table(out, file.path(output_dir, "monocle3_pseudotime.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cds, file.path(output_dir, "monocle3_cds.rds"))

write_status(
  "ok",
  "Monocle3 learn_graph + order_cells completed.",
  list(
    n_cells = ncol(cds),
    n_genes = nrow(cds),
    n_root_cells = length(root_cells),
    pseudotime_finite_cells = sum(is.finite(as.numeric(pt)))
  )
)
