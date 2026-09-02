#!/usr/bin/env Rscript

## Rebuild WT/KO network summaries from completed versioned RDS objects.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  equal_hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(equal_hit)) return(sub(paste0("^", flag, "="), "", equal_hit[1]))
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

rds_path <- get_arg("--rds", stop("--rds is required"))
output_path <- get_arg("--output", stop("--output is required"))
target_tf <- get_arg("--target-tf", stop("--target-tf is required"))
axis <- get_arg("--axis", "")
subset <- get_arg("--subset", "")
seed <- as.integer(get_arg("--seed", stop("--seed is required")))

edge_summary <- function(network, tf, label) {
  if (inherits(network, "Matrix") || inherits(network, "matrix")) {
    if (inherits(network, "matrix")) {
      values <- as.numeric(network)
      finite <- is.finite(values)
      nonzero <- finite & values != 0
      n_edges <- sum(nonzero)
      mean_abs <- if (n_edges) mean(abs(values[nonzero])) else 0
      max_abs <- if (n_edges) max(abs(values[nonzero])) else 0
    } else {
      sparse <- as(network, "dgCMatrix")
      edges <- Matrix::summary(sparse)
      n_edges <- nrow(edges)
      mean_abs <- if (n_edges) mean(abs(edges$x)) else 0
      max_abs <- if (n_edges) max(abs(edges$x)) else 0
    }
    row_names <- rownames(network)
    col_names <- colnames(network)
    n_nodes <- length(unique(c(row_names, col_names)))
    if (!n_nodes) n_nodes <- max(nrow(network), ncol(network))
    return(data.frame(tf = tf, axis = axis, subset = subset, seed = seed,
                      network = label, n_nodes = n_nodes, n_edges = n_edges,
                      mean_abs_weight = mean_abs, max_abs_weight = max_abs,
                      stringsAsFactors = FALSE))
  }
  data.frame()
}

recursive_summary <- function(item, tf, label) {
  if (is.list(item) && !is.data.frame(item) && !inherits(item, "Matrix")) {
    rows <- lapply(names(item), function(child) recursive_summary(item[[child]], tf, paste(label, child, sep = ".")))
    rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
    if (length(rows)) return(do.call(rbind, rows))
    return(data.frame())
  }
  edge_summary(item, tf, label)
}

result <- readRDS(rds_path)
if (is.null(result$tensorNetworks)) stop("tensorNetworks is missing from RDS")
rows <- lapply(names(result$tensorNetworks), function(name) recursive_summary(result$tensorNetworks[[name]], target_tf, paste0("tensorNetworks.", name)))
rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
summary <- if (length(rows)) do.call(rbind, rows) else data.frame(
  tf = target_tf, axis = axis, subset = subset, seed = seed,
  network = NA_character_, n_nodes = NA_integer_, n_edges = NA_integer_,
  mean_abs_weight = NA_real_, max_abs_weight = NA_real_
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(summary, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Network summary rebuilt:", output_path, "\n")
