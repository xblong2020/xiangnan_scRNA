#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

co <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "celloracle_module6_8_top_gene_delta_by_state.tsv.gz"))
sct <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "sctenifoldknk_module7_2_main_strict_perturbation_genes.tsv"))
ext <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "module8_tf_target_signature_genes.tsv"))
pathways <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "sctenifoldknk_module7_4_focus_tf_top_pathways.tsv"))
background <- sort(unique(as.character(sct$gene)))
focus <- c("HNF4A", "PPARA", "EGR1", "CEBPB", "AP-1 aggregate", "SOX4")

make_set <- function(tf_name) {
  members <- if (tf_name == "AP-1 aggregate") FIGURE6_AP1_MEMBERS else tf_name
  a <- co[tf %in% members & gene %in% background,
    .(score = max(abs_mean_delta_x, na.rm = TRUE)), by = gene][order(-score)][1:min(.N, 100)]
  b <- sct[tf %in% members & gene %in% background & p.adj < .05,
    .(score = max(distance, na.rm = TRUE), fdr = min(p.adj, na.rm = TRUE)), by = gene][order(fdr, -score)][1:min(.N, 100)]
  c <- ext[tf %in% members & gene %in% background][order(rank)][1:min(.N, 50)]
  rbindlist(list(
    a[, .(gene, evidence_source = "CellOracle top predicted changes", evidence_value = score)],
    b[, .(gene, evidence_source = "scTenifoldKnk FDR<0.05", evidence_value = score)],
    c[, .(gene, evidence_source = "Frozen external/signature target", evidence_value = as.numeric(rank))]
  ), fill = TRUE)[, .(evidence_sources = paste(sort(unique(evidence_source)), collapse = ";"),
    n_evidence_sources = uniqueN(evidence_source)), by = gene][, perturbation := tf_name]
}
sets <- rbindlist(lapply(focus, make_set), fill = TRUE)
sets[, `:=`(background_n = length(background), background_definition = "module7.2 main_strict 3000-gene network")]
figure6_fwrite(sets, file.path(FIGURE6_METADATA_DIR, "figure6e_gene_sets.tsv"), compress = FALSE)

pairs <- CJ(set_1 = focus, set_2 = focus)
pairs[, jaccard := mapply(function(a, b) figure6_jaccard(sets[perturbation == a, gene], sets[perturbation == b, gene]), set_1, set_2)]
pairs[, overlap_n := mapply(function(a, b) length(intersect(sets[perturbation == a, gene], sets[perturbation == b, gene])), set_1, set_2)]
pairs[, set_type := "gene"]
figure6_fwrite(pairs, file.path(FIGURE6_METADATA_DIR, "figure6e_jaccard_matrix.tsv"), compress = FALSE)

sig_pathways <- pathways[is.finite(p.adjust) & p.adjust < .05]
path_pairs <- CJ(set_1 = focus, set_2 = focus)
path_pairs[, `:=`(
  jaccard = mapply(function(a, b) figure6_jaccard(
    sig_pathways[tf %in% if (a == "AP-1 aggregate") FIGURE6_AP1_MEMBERS else a, term_id],
    sig_pathways[tf %in% if (b == "AP-1 aggregate") FIGURE6_AP1_MEMBERS else b, term_id]), set_1, set_2),
  n_fdr_pathways_set_1 = vapply(set_1, function(a) uniqueN(sig_pathways[tf %in% if (a == "AP-1 aggregate") FIGURE6_AP1_MEMBERS else a, term_id]), integer(1)),
  n_fdr_pathways_set_2 = vapply(set_2, function(a) uniqueN(sig_pathways[tf %in% if (a == "AP-1 aggregate") FIGURE6_AP1_MEMBERS else a, term_id]), integer(1))
)]
path_pairs[, status := fifelse(n_fdr_pathways_set_1 == 0 | n_fdr_pathways_set_2 == 0, "not_estimable_no_FDR_pathways", "estimable")]
figure6_fwrite(path_pairs, file.path(FIGURE6_METADATA_DIR, "figure6e_pathway_similarity.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6E", background = "module7.2 main_strict network", background_n = length(background),
  perturbations = as.list(focus), gene_set_rule = "union of frozen CellOracle top-100, scTenifoldKnk FDR<0.05 top-100 and frozen target-signature top-50",
  n_fdr_pathways = nrow(sig_pathways), pathway_result = if (nrow(sig_pathways)) "available" else "no FDR-significant pathway; similarity not estimable",
  guardrail = "No pathway category was prespecified or inferred from nonsignificant enrichment."
), file.path(FIGURE6_METADATA_DIR, "figure6e_overlap_report.json"))

