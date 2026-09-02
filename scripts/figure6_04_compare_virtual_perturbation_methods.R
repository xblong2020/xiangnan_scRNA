#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

co <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "celloracle_module6_8_top_gene_delta_by_state.tsv.gz"))
sct <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "sctenifoldknk_module7_2_main_strict_perturbation_genes.tsv"))
tfs <- intersect(unique(co$tf), unique(sct$tf))
rows <- lapply(seq_along(tfs), function(i) {
  tf_i <- tfs[i]
  a <- co[tf == tf_i, .(celloracle_magnitude = mean(abs_mean_delta_x, na.rm = TRUE)), by = gene]
  b <- sct[tf == tf_i, .(sct_magnitude = max(distance, na.rm = TRUE), sct_fdr = min(p.adj, na.rm = TRUE)), by = gene]
  a[, co_rank := frank(-celloracle_magnitude, ties.method = "average")]
  b[, sct_rank := frank(-sct_magnitude, ties.method = "average")]
  m <- merge(a, b, by = "gene")
  top_co <- a[order(co_rank)][1:min(50, .N), gene]
  top_sct <- b[order(sct_rank)][1:min(50, .N), gene]
  rho <- if (nrow(m) >= 5) suppressWarnings(cor(m$co_rank, m$sct_rank, method = "spearman")) else NA_real_
  data.table(
    tf = tf_i, axis = figure6_axis_for_tf(tf_i), n_shared_ranked_genes = nrow(m), spearman_rank_correlation = rho,
    top50_jaccard = figure6_jaccard(top_co, top_sct), top50_overlap = length(intersect(top_co, top_sct)),
    sct_fdr_gene_overlap = length(intersect(top_co, b[sct_fdr < .05, gene])),
    sign_concordance = NA_real_, direction_interpretability = "magnitude/rank only",
    dataset_stability = NA_real_, signed_celloracle = TRUE, signed_sctenifold = FALSE
  )
})
out <- rbindlist(rows)
figure6_fwrite(out, file.path(FIGURE6_METADATA_DIR, "figure6d_cross_method_concordance.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6D", comparison = "CellOracle absolute predicted-expression change versus scTenifoldKnk manifold distance",
  scTenifoldKnk_directionality = "unsigned magnitude only", top_k = 50, n_tfs = nrow(out),
  prohibited_interpretation = "Z and FC columns are not used as signed perturbation effects"
), file.path(FIGURE6_METADATA_DIR, "figure6d_cross_method_report.json"))

