#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_paths(root)

target_path <- file.path(root, "metadata/driver/module8_tf_target_signature_genes.tsv")
celloracle_path <- file.path(root, "metadata/driver/celloracle_module6_7_grn_links_filtered.tsv.gz")
targets <- fread(target_path)
links <- fread(celloracle_path, showProgress = FALSE)

axis_tfs <- list(
  identity_loss = c("HNF4A", "PPARA"),
  stress_transition = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "ATF3", "CEBPB", "EGR1"),
  sox4_stabilization = "SOX4"
)
directions <- c(identity_loss = "decreasing_identity_retention", stress_transition = "increasing", sox4_stabilization = "increasing")

core <- rbindlist(lapply(names(axis_tfs), function(axis) data.table(
  axis = axis, gene = axis_tfs[[axis]], source = "predefined_core_TF", source_figure = "Figures 2-4",
  selection_rule = "predefined axis TF; no pseudotime selection", direction = directions[[axis]],
  primary_or_sensitivity = "primary", version = "figure5_contract_v1", file_path = "Figure 5 task specification",
  rank = seq_along(axis_tfs[[axis]]), entity_type = "TF expression"
)))

sct <- rbindlist(lapply(names(axis_tfs), function(axis) {
  axis_name <- axis
  sub <- targets[tf %chin% axis_tfs[[axis_name]] & signature_class == "main"]
  if (!nrow(sub)) return(NULL)
  sub[, .(axis = axis_name, gene = toupper(gene), source = "scTenifoldKnk_frozen_target", source_figure = "Figures 2-4 / module8",
          selection_rule = paste0("frozen module8 target rank; TF in ", paste(axis_tfs[[axis_name]], collapse = "/")), direction = directions[[axis_name]],
          primary_or_sensitivity = "sensitivity", version = unname(tools::md5sum(target_path)), file_path = normalizePath(target_path, winslash = "/"),
          rank = as.integer(rank), entity_type = "target gene")]
}))

links[, source := toupper(source)]
links[, target := toupper(target)]
co <- rbindlist(lapply(names(axis_tfs), function(axis) {
  axis_name <- axis
  sub <- links[source %chin% axis_tfs[[axis_name]] & is.finite(coef_abs)]
  if (!nrow(sub)) return(NULL)
  ranked <- sub[, .(coef_abs = max(coef_abs, na.rm = TRUE), min_p = min(p, na.rm = TRUE)), by = .(source, target)][order(-coef_abs)]
  ranked <- ranked[, head(.SD, 10L), by = source]
  ranked[, .(axis = axis_name, gene = target, source = "CellOracle_filtered_link", source_figure = "Figures 2-4 / module6",
             selection_rule = "top 10 frozen filtered targets per TF by absolute coefficient; no pseudotime selection", direction = directions[[axis_name]],
             primary_or_sensitivity = "sensitivity", version = unname(tools::md5sum(celloracle_path)), file_path = normalizePath(celloracle_path, winslash = "/"),
             rank = frank(-coef_abs, ties.method = "first"), entity_type = "target gene")]
}))

intersection <- merge(unique(sct[, .(axis, gene)]), unique(co[, .(axis, gene)]), by = c("axis", "gene"))
if (nrow(intersection)) {
  intersection[, rank := seq_len(.N), by = axis]
  intersection[, `:=`(
    source = "CellOracle_scTenifold_intersection", source_figure = "Figures 2-4",
    selection_rule = "present in both frozen CellOracle and scTenifoldKnk target sets", direction = directions[axis],
    primary_or_sensitivity = "sensitivity", version = "intersection_v1",
    file_path = paste(normalizePath(target_path, winslash = "/"), normalizePath(celloracle_path, winslash = "/"), sep = ";"),
    entity_type = "target gene"
  )]
}

audit <- unique(rbindlist(list(core, sct, co, intersection), fill = TRUE), by = c("axis", "gene", "source"))
setorder(audit, axis, primary_or_sensitivity, source, rank, gene)
figure5_write_tsv(audit, file.path(paths$metadata, "figure5_frozen_signature_audit.tsv"))

# Heatmap selection is frozen by source rank and availability is checked later against the H5AD gene universe.
heatmap_targets <- rbindlist(lapply(names(axis_tfs), function(axis) {
  axis_name <- axis
  preferred <- audit[axis == axis_name & source == "CellOracle_scTenifold_intersection"][order(rank)]
  if (nrow(preferred) < 6L) preferred <- unique(rbindlist(list(preferred, audit[axis == axis_name & source == "scTenifoldKnk_frozen_target"][order(rank)]), fill = TRUE), by = "gene")
  head(preferred, 6L)
}))
heatmap_manifest <- unique(rbindlist(list(core, heatmap_targets), fill = TRUE), by = c("axis", "gene", "entity_type"))
figure5_write_tsv(heatmap_manifest, file.path(paths$metadata, "figure5_heatmap_gene_manifest_preavailability.tsv"))

generic_stress <- c("JUN", "JUNB", "JUND", "FOS", "FOSB", "ATF3", "HSPA1A", "HSPA1B", "HSPA6", "DDIT3")
cell_cycle <- c("MKI67", "TOP2A", "BIRC5", "UBE2C", "CENPF", "TYMS", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7")
gene_sets <- list(
  tf_identity_retention = axis_tfs$identity_loss,
  tf_stress_transition = axis_tfs$stress_transition,
  tf_sox4_stabilization = axis_tfs$sox4_stabilization,
  celloracle_identity_retention = unique(co[axis == "identity_loss"]$gene),
  celloracle_stress_transition = unique(co[axis == "stress_transition"]$gene),
  celloracle_sox4_stabilization = unique(co[axis == "sox4_stabilization"]$gene),
  intersection_identity_retention = unique(intersection[axis == "identity_loss"]$gene),
  intersection_stress_transition = unique(intersection[axis == "stress_transition"]$gene),
  intersection_sox4_stabilization = unique(intersection[axis == "sox4_stabilization"]$gene),
  identity_no_cell_cycle = setdiff(unique(c(axis_tfs$identity_loss, sct[axis == "identity_loss"]$gene)), cell_cycle),
  stress_no_cell_cycle = setdiff(unique(c(axis_tfs$stress_transition, sct[axis == "stress_transition"]$gene)), cell_cycle),
  sox4_no_cell_cycle = setdiff(unique(c(axis_tfs$sox4_stabilization, sct[axis == "sox4_stabilization"]$gene)), cell_cycle),
  stress_no_generic = setdiff(unique(c(axis_tfs$stress_transition, sct[axis == "stress_transition"]$gene)), generic_stress)
)
saveRDS(gene_sets, file.path(paths$processed, "figure5_frozen_gene_sets.rds"))
axis_counts <- audit[, .(n_genes = uniqueN(gene)), by = axis]
figure5_write_json(list(n_genes_by_axis = as.list(setNames(axis_counts$n_genes, axis_counts$axis)),
                        source_files = c(target_path, celloracle_path), selection_used_pseudotime = FALSE),
                   file.path(paths$metadata, "figure5_frozen_signature_audit.json"))
message("Frozen signature audit complete: ", nrow(audit), " source rows.")
