#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_brd_core <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  keep <- !is.na(x) & grepl("BRD-[A-Z][0-9]{8}", x, perl = TRUE)
  out[keep] <- sub("^.*?(BRD-[A-Z][0-9]{8}).*$", "\\1", x[keep], perl = TRUE)
  out
}

figure8_v2_network_score <- function(direct_target, one_step, pathway_overlap, axis_compatibility) {
  values <- c(direct_target, one_step, pathway_overlap, axis_compatibility)
  values[!is.finite(values)] <- 0
  pmax(0, pmin(1, sum(values * c(0.40, 0.25, 0.20, 0.15))))
}

figure8_v2_axis_compatibility <- function(target, moa) {
  target <- toupper(trimws(as.character(target)))
  moa_lower <- tolower(as.character(moa))
  positive_action <- grepl("agonist|activator|stabiliz|induc|enhanc", moa_lower)
  negative_action <- grepl("inhibitor|antagonist|blocker|suppress|degrad", moa_lower)
  if (target %in% c("HNF4A", "PPARA", "HLF")) {
    return(list(axis = "axis_A_identity", score = as.numeric(positive_action), compatibility = if (positive_action) "compatible" else if (negative_action) "opposing" else "direction_unknown"))
  }
  if (target %in% c("JUN", "JUNB", "JUND", "FOS", "CEBPB", "EGR1", "ATF3")) {
    return(list(axis = "axis_B_stress", score = as.numeric(negative_action), compatibility = if (negative_action) "compatible" else if (positive_action) "opposing" else "direction_unknown"))
  }
  if (target == "SOX4") {
    return(list(axis = "axis_C_sox4", score = as.numeric(negative_action), compatibility = if (negative_action) "compatible" else if (positive_action) "opposing" else "direction_unknown"))
  }
  list(axis = NA_character_, score = NA_real_, compatibility = "not_axis_target")
}

figure8_v2_split_targets <- function(x) {
  x <- as.character(x %||% "")
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  targets <- unlist(strsplit(x, "[,;|/]", perl = TRUE), use.names = FALSE)
  targets <- toupper(trimws(targets))
  sort(unique(targets[nzchar(targets) & targets != "NA"]))
}

figure8_v2_target_edges <- function(annotation) {
  x <- as.data.table(copy(annotation))[!is.na(curated_targets) & curated_targets != ""]
  if (!nrow(x)) return(data.table(candidate = character(), BRD_ID = character(), target = character(), edge_type = character(), source = character(), confidence = character()))
  x[, row_id__ := .I]
  edges <- x[, {
    targets <- figure8_v2_split_targets(curated_targets[[1]])
    if (!length(targets)) NULL else data.table(
      candidate = canonical_name[[1]], BRD_ID = BRD_ID[[1]], target = targets,
      edge_type = "curated_target", source = source[[1]], confidence = annotation_confidence[[1]]
    )
  }, by = row_id__]
  edges[, row_id__ := NULL]
  edges
}

figure8_v2_moa_network <- function() {
  raw_dir <- file.path(FIGURE8_V2_DATA, "external_raw")
  primary <- fread(file.path(raw_dir, "Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv"))
  primary_ann <- primary[, .(
    brd_core = figure8_v2_brd_core(IDs),
    source_name = as.character(Drug.Name),
    curated_MoA = as.character(MOA),
    curated_targets = as.character(repurposing_target),
    clinical_status = NA_character_,
    approved_status = NA_character_,
    synonyms = as.character(Synonyms),
    annotation_source = "Broad PRISM Repurposing Public 23Q2 v4 compound list",
    annotation_source_id = as.character(IDs),
    annotation_confidence = "high_exact_BRD"
  )]

  secondary <- fread(file.path(raw_dir, "secondary-screen-replicate-collapsed-treatment-info.csv"))
  secondary_ann <- unique(secondary[, .(
    brd_core = figure8_v2_brd_core(broad_id),
    source_name = as.character(name),
    curated_MoA = as.character(moa),
    curated_targets = as.character(target),
    clinical_status = as.character(phase),
    approved_status = fifelse(grepl("launched|approved", phase, ignore.case = TRUE), "approved_or_launched", "not_confirmed_approved"),
    synonyms = NA_character_,
    annotation_source = "Broad PRISM Repurposing 19Q4 v4 secondary treatment metadata",
    annotation_source_id = as.character(broad_id),
    annotation_confidence = "high_exact_BRD"
  )])
  annotations <- rbindlist(list(primary_ann, secondary_ann), fill = TRUE)
  annotations <- annotations[!is.na(brd_core)]
  collapse_values <- function(x) {
    values <- sort(unique(trimws(as.character(x[!is.na(x) & as.character(x) != ""]))))
    if (length(values)) paste(values, collapse = ";") else NA_character_
  }
  curated <- annotations[, .(
    source_name = collapse_values(source_name),
    curated_MoA = collapse_values(curated_MoA),
    curated_targets = collapse_values(curated_targets),
    clinical_status = collapse_values(clinical_status),
    approved_status = collapse_values(approved_status),
    synonyms = collapse_values(synonyms),
    source = collapse_values(annotation_source),
    source_identifiers = collapse_values(annotation_source_id),
    annotation_confidence = "high_exact_BRD",
    mapping_conflict = uniqueN(tolower(source_name[!is.na(source_name)])) > 1L | uniqueN(tolower(curated_MoA[!is.na(curated_MoA)])) > 1L
  ), by = brd_core]

  ranking <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_full_ranking.tsv.gz"))
  ranking[, brd_core := figure8_v2_brd_core(compound)]
  annotation <- merge(ranking[, .(
    BRD_ID = compound, brd_core, canonical_name, InChIKey = inchi_key,
    PubChem_CID = pubchem_cid, v2_primary_rank, candidate_analysis_universe
  )], curated, by = "brd_core", all.x = TRUE)
  annotation[, `:=`(
    curated_annotation_available = !is.na(curated_MoA) | !is.na(curated_targets),
    target_class = NA_character_,
    retrieval_date = "2026-08-18",
    inferred_mechanism_used_as_curated = FALSE
  )]
  figure8_v2_write_tsv(annotation, "figure8_v2_compound_moa_target_annotation.tsv")

  core_tfs <- c("HNF4A", "PPARA", "HLF", "JUN", "JUNB", "JUND", "FOS", "CEBPB", "EGR1", "ATF3", "SOX4")
  one_step_genes <- unique(fread(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure6_directional_network/figure6e_gene_sets.tsv"))$gene)
  one_step_genes <- toupper(as.character(one_step_genes))
  pathway_genes <- unique(fread(file.path(FIGURE8_V2_ROOT, "metadata/driver/module8_pathway_signature_genes.tsv"))$gene)
  pathway_genes <- toupper(as.character(pathway_genes))

  network_rows <- lapply(seq_len(nrow(annotation)), function(idx) {
    row <- annotation[idx]
    targets <- figure8_v2_split_targets(row$curated_targets)
    if (!length(targets)) {
      return(data.table(
        BRD_ID = row$BRD_ID, direct_target_score = NA_real_, one_step_network_proximity = NA_real_,
        pathway_overlap = NA_real_, axis_direction_compatibility = NA_real_, network_consistency_score = NA_real_,
        network_evidence_coverage = 0, direct_network_targets = NA_character_, one_step_targets = NA_character_,
        pathway_targets = NA_character_, compatible_axes = NA_character_,
        network_consistent_inferred_mechanism = NA_character_
      ))
    }
    direct <- targets %in% core_tfs
    one_step <- targets %in% one_step_genes
    pathway <- targets %in% pathway_genes
    compatibility <- lapply(targets, figure8_v2_axis_compatibility, moa = row$curated_MoA %||% "")
    compatibility_scores <- vapply(compatibility, `[[`, numeric(1), "score")
    compatible_axes <- unique(vapply(compatibility[is.finite(compatibility_scores) & compatibility_scores > 0], `[[`, character(1), "axis"))
    components <- c(
      direct_target = mean(direct), one_step = mean(one_step), pathway_overlap = mean(pathway),
      axis_compatibility = if (any(is.finite(compatibility_scores))) max(compatibility_scores, na.rm = TRUE) else NA_real_
    )
    data.table(
      BRD_ID = row$BRD_ID,
      direct_target_score = components[["direct_target"]],
      one_step_network_proximity = components[["one_step"]],
      pathway_overlap = components[["pathway_overlap"]],
      axis_direction_compatibility = components[["axis_compatibility"]],
      network_consistency_score = figure8_v2_network_score(components[[1]], components[[2]], components[[3]], components[[4]]),
      network_evidence_coverage = mean(is.finite(components)),
      direct_network_targets = if (any(direct)) paste(targets[direct], collapse = ";") else NA_character_,
      one_step_targets = if (any(one_step)) paste(targets[one_step], collapse = ";") else NA_character_,
      pathway_targets = if (any(pathway)) paste(targets[pathway], collapse = ";") else NA_character_,
      compatible_axes = if (length(compatible_axes)) paste(sort(compatible_axes), collapse = ";") else NA_character_,
      network_consistent_inferred_mechanism = if (length(compatible_axes)) paste0("network-consistent with ", paste(sort(compatible_axes), collapse = ", ")) else "no direction-specific network compatibility"
    )
  })
  network <- rbindlist(network_rows, fill = TRUE)
  network <- merge(network, annotation[, .(BRD_ID, canonical_name, curated_MoA, curated_targets, annotation_confidence, mapping_conflict)], by = "BRD_ID", all.x = TRUE)
  figure8_v2_write_tsv(network, "figure8_v2_network_consistency.tsv")

  edge_rows <- figure8_v2_target_edges(annotation)
  figure8_v2_write_tsv(edge_rows, "figure8_v2_candidate_target_edges.tsv")

  figure8_v2_write_json(list(
    module = "figure8_v2_moa_network", status = "completed",
    n_compounds = nrow(annotation),
    moa_annotation_rate = mean(!is.na(annotation$curated_MoA)),
    target_annotation_rate = mean(!is.na(annotation$curated_targets)),
    network_consistent_count = sum(network$network_consistency_score >= 0.50, na.rm = TRUE),
    score_definition = "0.40 direct target + 0.25 one-step proximity + 0.20 pathway overlap + 0.15 axis-direction compatibility",
    provenance_boundary = "Curated MoA/targets come only from official compound metadata; Figure 6 evidence is separately labelled inferred network consistency."
  ), "figure8_v2_moa_network_report.json")
  invisible(list(annotation = annotation, network = network))
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_moa_network()
  cat("FIGURE8_V2_MOA annotated=", sum(result$annotation$curated_annotation_available), "/", nrow(result$annotation), " network_consistent=", sum(result$network$network_consistency_score >= 0.50, na.rm = TRUE), "\n", sep = "")
}
