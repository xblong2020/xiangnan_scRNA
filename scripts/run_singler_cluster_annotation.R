suppressPackageStartupMessages({
  library(SingleR)
  library(celldex)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
flags <- args[startsWith(args, "--")]
positional <- args[!startsWith(args, "--")]
reuse_existing <- "--reuse-existing" %in% flags
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (length(positional) >= 1) {
  matrix_path <- positional[[1]]
} else {
  matrix_path <- file.path(root, "metadata/celltype/singler_cluster_logcpm.tsv.gz")
}
if (length(positional) >= 2) {
  out_dir <- positional[[2]]
} else {
  out_dir <- file.path(root, "metadata/celltype")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  hpca_main = file.path(out_dir, "singler_hpca_main_by_cluster.tsv"),
  hpca_fine = file.path(out_dir, "singler_hpca_fine_by_cluster.tsv"),
  blueprint_main = file.path(out_dir, "singler_blueprint_main_by_cluster.tsv"),
  blueprint_fine = file.path(out_dir, "singler_blueprint_fine_by_cluster.tsv")
)

if (!reuse_existing || !all(file.exists(unlist(paths)))) {
  message("READ ", matrix_path)
  expr <- read.delim(matrix_path, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  message("MATRIX genes=", nrow(expr), " clusters=", ncol(expr))
} else {
  expr <- NULL
  message("REUSE existing SingleR per-reference outputs")
}

add_top_score <- function(out, pred = NULL) {
  if (!("singler_top_score" %in% colnames(out))) {
    if (!is.null(pred)) {
      score_matrix <- as.matrix(pred$scores)
      out$singler_top_score <- apply(score_matrix, 1, max, na.rm = TRUE)
    } else {
      score_cols <- grep("^scores\\.", colnames(out), value = TRUE)
      if (length(score_cols) == 0) {
        out$singler_top_score <- NA_real_
      } else {
        out$singler_top_score <- apply(out[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
      }
    }
  }
  out
}

run_one <- function(ref, ref_name, labels) {
  common <- intersect(rownames(expr), rownames(ref))
  if (length(common) < 1000) {
    stop("Too few common genes for ", ref_name, ": ", length(common))
  }
  message("RUN ", ref_name, " common_genes=", length(common))
  pred <- SingleR(test = expr[common, , drop = FALSE], ref = ref[common, ], labels = labels)
  out <- as.data.frame(pred)
  out <- add_top_score(out, pred)
  out$singler_column <- rownames(out)
  out$reference <- ref_name
  out$common_genes <- length(common)
  out
}

if (reuse_existing && all(file.exists(unlist(paths)))) {
  hpca_main <- add_top_score(read.delim(paths$hpca_main, check.names = FALSE))
  hpca_fine <- add_top_score(read.delim(paths$hpca_fine, check.names = FALSE))
  bp_main <- add_top_score(read.delim(paths$blueprint_main, check.names = FALSE))
  bp_fine <- add_top_score(read.delim(paths$blueprint_fine, check.names = FALSE))
} else {
  hpca <- HumanPrimaryCellAtlasData()
  bp <- BlueprintEncodeData()

  hpca_main <- run_one(hpca, "HumanPrimaryCellAtlasData_label.main", hpca$label.main)
  hpca_fine <- run_one(hpca, "HumanPrimaryCellAtlasData_label.fine", hpca$label.fine)
  bp_main <- run_one(bp, "BlueprintEncodeData_label.main", bp$label.main)
  bp_fine <- run_one(bp, "BlueprintEncodeData_label.fine", bp$label.fine)

  write.table(hpca_main, paths$hpca_main, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(hpca_fine, paths$hpca_fine, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bp_main, paths$blueprint_main, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bp_fine, paths$blueprint_fine, sep = "\t", quote = FALSE, row.names = FALSE)
}

select_pred <- function(x, prefix) {
  out <- x[, c("singler_column", "labels", "pruned.labels", "singler_top_score", "delta.next", "common_genes")]
  colnames(out) <- c(
    "singler_column",
    paste0(prefix, "_label"),
    paste0(prefix, "_pruned_label"),
    paste0(prefix, "_top_score"),
    paste0(prefix, "_delta_next"),
    paste0(prefix, "_common_genes")
  )
  out
}

singler_summary <- Reduce(
  function(left, right) merge(left, right, by = "singler_column", all = TRUE, sort = FALSE),
  list(
    select_pred(hpca_main, "singler_hpca_main"),
    select_pred(hpca_fine, "singler_hpca_fine"),
    select_pred(bp_main, "singler_blueprint_main"),
    select_pred(bp_fine, "singler_blueprint_fine")
  )
)

celltypist_path <- file.path(out_dir, "celltypist_major_by_leiden.tsv")
cluster_meta_path <- file.path(out_dir, "singler_cluster_metadata.tsv")
if (file.exists(celltypist_path)) {
  cluster_meta <- read.delim(celltypist_path, check.names = FALSE)
  cluster_meta$singler_column <- paste0("cluster_", cluster_meta$leiden_scvi)
  if ("major_celltype" %in% colnames(cluster_meta)) {
    colnames(cluster_meta)[colnames(cluster_meta) == "major_celltype"] <- "celltypist_major"
  }
  if ("major_celltype_fraction" %in% colnames(cluster_meta)) {
    colnames(cluster_meta)[colnames(cluster_meta) == "major_celltype_fraction"] <- "celltypist_major_fraction"
  }
  merged <- merge(cluster_meta, singler_summary, by = "singler_column", all.y = TRUE, sort = FALSE)
} else if (file.exists(cluster_meta_path)) {
  cluster_meta <- read.delim(cluster_meta_path, check.names = FALSE)
  merged <- merge(cluster_meta, singler_summary, by = "singler_column", all.y = TRUE, sort = FALSE)
} else {
  merged <- singler_summary
}
write.table(merged, file.path(out_dir, "singler_combined_by_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

report <- list(
  input = normalizePath(matrix_path, winslash = "/", mustWork = FALSE),
  output = normalizePath(file.path(out_dir, "singler_combined_by_cluster.tsv"), winslash = "/", mustWork = FALSE),
  reused_existing_reference_outputs = reuse_existing && all(file.exists(unlist(paths))),
  n_genes = if (is.null(expr)) jsonlite::read_json(file.path(out_dir, "singler_cluster_expression_report.json"))$n_genes else nrow(expr),
  n_clusters = nrow(merged),
  hpca_common_genes = unique(hpca_main$common_genes),
  blueprint_common_genes = unique(bp_main$common_genes),
  SingleR_version = as.character(packageVersion("SingleR")),
  celldex_version = as.character(packageVersion("celldex"))
)
writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "singler_cluster_report.json"))
message("WROTE ", file.path(out_dir, "singler_combined_by_cluster.tsv"))
