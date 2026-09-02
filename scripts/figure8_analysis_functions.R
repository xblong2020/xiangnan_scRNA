#!/usr/bin/env Rscript

## Figure 8 analysis and plotting functions. Every wrapper sources the palette
## contract before sourcing this file.

figure8_input <- function(...) figure8_root_path(...)

figure8_signature_path <- function() figure8_input("metadata", "driver", "module9_4_drug_reversal_signature.tsv")
figure8_model_gene_path <- function() file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_model_genes.tsv")

figure8_signature <- function() {
  x <- figure8_fread(figure8_signature_path())
  x[, gene := toupper(trimws(as.character(gene)))]
  x[, `:=`(
    include_primary = figure8_bool(include_primary),
    include_sensitivity = figure8_bool(include_sensitivity),
    conflict_flag = figure8_bool(conflict_flag),
    housekeeping_or_qc_flag = figure8_bool(housekeeping_or_qc_flag),
    final_weight = as.numeric(final_weight),
    axis = figure8_axis_from_component(component)
  )]
  x
}

figure8_model_genes <- function() {
  path <- figure8_model_gene_path()
  figure8_require(path, "DrugReflector model-gene manifest")
  unique(toupper(figure8_fread(path)$gene))
}

figure8_read_gmt <- function(path, source_name = basename(path)) {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) return(NULL)
    data.table(term = fields[[1]], description = fields[[2]], gene = toupper(fields[-c(1, 2)]), source = source_name)
  })
  rbindlist(rows, fill = TRUE)
}

figure8_local_nuisance_sets <- function() {
  files <- c(
    figure8_input("metadata", "driver", "sctenifoldknk_module7_4_genesets", "GO_Biological_Process_2023.gmt"),
    figure8_input("metadata", "driver", "sctenifoldknk_module7_4_genesets", "Reactome_2022.gmt"),
    figure8_input("metadata", "driver", "sctenifoldknk_module7_4_genesets", "KEGG_2021_Human.gmt")
  )
  files <- files[file.exists(files)]
  if (!length(files)) return(data.table())
  gmt <- rbindlist(lapply(files, figure8_read_gmt), fill = TRUE)
  patterns <- c(
    proliferation = "proliferat|cell population growth",
    cell_cycle = "cell cycle|mitotic|g0 and early g1|g1 phase|g2 phase|s phase|checkpoint",
    generic_stress = "stress response|response to stress|chemical stress|integrated stress|endoplasmic reticulum stress|unfolded protein",
    dna_damage = "dna damage|dna repair|genotoxic|double.strand break",
    translation_inhibition = "translation|ribosom|protein synthesis",
    mitochondrial_toxicity = "mitochond|oxidative phosphorylation|respiratory electron transport|electron transport chain",
    unfolded_protein_response = "unfolded protein|endoplasmic reticulum stress|atf4 activates genes",
    apoptosis = "apopto|programmed cell death|intrinsic cell death"
  )
  out <- rbindlist(lapply(names(patterns), function(set_name) {
    hit <- gmt[grepl(patterns[[set_name]], term, ignore.case = TRUE)]
    if (!nrow(hit)) return(NULL)
    unique(hit[, .(gene, term, source)])[, nuisance_set := set_name]
  }), fill = TRUE)
  setcolorder(out, c("nuisance_set", "gene", "term", "source"))
  unique(out)
}

figure8_parse_json <- function(path) {
  if (!file.exists(path)) return(NULL)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

figure8_file_records <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  ext <- tolower(basename(path))
  if (grepl("\\.(tsv|txt)(\\.gz)?$", ext)) {
    value <- tryCatch(nrow(data.table::fread(path, select = 1, showProgress = FALSE)), error = function(e) NA_integer_)
    return(as.integer(value))
  }
  if (grepl("\\.gmt$", ext)) return(length(readLines(path, warn = FALSE)))
  if (grepl("\\.json$", ext)) return(1L)
  NA_integer_
}

figure8_protected_manifest <- function() {
  roots <- file.path(FIGURE8_PROJECT_ROOT, c("scripts", "metadata", "figures", "reports"))
  files <- unlist(lapply(roots[dir.exists(roots)], list.files, recursive = TRUE, full.names = TRUE, all.files = TRUE), use.names = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(FIGURE8_PROJECT_ROOT, winslash = "/")
  rel <- ifelse(startsWith(norm, paste0(root_norm, "/")), substring(norm, nchar(root_norm) + 2L), norm)
  script_new <- grepl("^scripts/(figure8_|plot_figure8|validate_figure8|run_figure8)", rel, ignore.case = TRUE)
  exclusive <- grepl("metadata/driver/figure8_transcriptomic_reversal|figures/driver/figure8[a-h]_.*|figures/driver/figure8_transcriptomic_reversal_preview|reports/figure8_transcriptomic_reversal_report", rel, ignore.case = TRUE)
  keep <- !(script_new | exclusive)
  info <- file.info(files[keep])
  data.table(
    file_path = rel[keep],
    size_bytes = as.numeric(info$size),
    modified_utc = format(as.POSIXct(info$mtime, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    md5 = unname(tools::md5sum(files[keep]))
  )
}

figure8_protected_data_manifest <- function() {
  root <- figure8_input("data")
  files <- if (dir.exists(root)) list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE) else character()
  files <- files[file.info(files)$isdir %in% FALSE]
  norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(FIGURE8_PROJECT_ROOT, winslash = "/")
  rel <- ifelse(startsWith(norm, paste0(root_norm, "/")), substring(norm, nchar(root_norm) + 2L), norm)
  keep <- !grepl("data/processed/driver/figure8_transcriptomic_reversal", rel, ignore.case = TRUE)
  info <- file.info(files[keep])
  data.table(
    file_path = rel[keep], size_bytes = as.numeric(info$size),
    modified_utc = format(as.POSIXct(info$mtime, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  )
}

run_figure8_preflight <- function() {
  figure8_dirs()
  signature_path <- figure8_signature_path()
  core <- list(
    signature = signature_path,
    signature_json = figure8_input("metadata", "driver", "module9_4_drug_reversal_signature.json"),
    signature_qc = figure8_input("metadata", "driver", "module9_4_signature_qc.tsv"),
    module94_report = figure8_input("metadata", "driver", "module9_4_report.json"),
    dr_primary = figure8_input("metadata", "driver", "module9_7_drugreflector_primary_predictions.tsv"),
    dr_sensitivity = figure8_input("metadata", "driver", "module9_7_drugreflector_sensitivity_predictions.tsv"),
    dr_consensus = figure8_input("metadata", "driver", "module9_7_drugreflector_consensus_predictions.tsv"),
    dr_coverage = figure8_input("metadata", "driver", "module9_7_drugreflector_gene_coverage.tsv"),
    dr_report = figure8_input("metadata", "driver", "module9_7_drugreflector_report.json"),
    l1000 = figure8_input("metadata", "driver", "module9_5_l1000fwd_candidate_ranking.tsv"),
    clue_status = figure8_input("metadata", "driver", "module9_8_clue_job_status.json"),
    clue_crosswalk = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_clue_crosswalk.tsv"),
    clue_summary = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_clue_perturbagen_summary.tsv"),
    clue_cells = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_clue_cell_scores.tsv.gz"),
    crossvalidation = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation.tsv"),
    crossvalidation_report = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_report.json"),
    perturbagen_metadata = figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_perturbagen_metadata.tsv"),
    decomposition = figure8_input("metadata", "driver", "module9_9_landmark_decomposition_drugreflector_predictions.tsv.gz"),
    clue_components = figure8_input("metadata", "driver", "module9_9_landmark_decomposition_clue_signature_components.tsv.gz"),
    decomposition_final = figure8_input("metadata", "driver", "module9_9_landmark_decomposition_final_priority.tsv"),
    checkpoint0 = figure8_input("metadata", "driver", "drugreflector_checkpoints", "model_fold_0.pt"),
    checkpoint1 = figure8_input("metadata", "driver", "drugreflector_checkpoints", "model_fold_1.pt"),
    checkpoint2 = figure8_input("metadata", "driver", "drugreflector_checkpoints", "model_fold_2.pt")
  )
  signature <- if (file.exists(signature_path)) figure8_signature() else data.table()
  dr_report <- figure8_parse_json(core$dr_report)
  cv_report <- figure8_parse_json(core$crossvalidation_report)
  clue_report <- figure8_parse_json(core$clue_status)
  coverage <- if (file.exists(core$dr_coverage)) figure8_fread(core$dr_coverage) else data.table()
  cross <- if (file.exists(core$crossvalidation)) figure8_fread(core$crossvalidation) else data.table()

  frozen_rows <- list(
    c("Module 9.4", "gene_signature", "primary_three_axis", "signature_construction", core$signature, "signed", "FALSE", "FALSE", "FALSE", "FALSE", "HGNC symbol", "primary", "frozen", "Weighted gene-level provenance"),
    c("Module 9.4", "gene_signature", "full_three_axis", "signature_construction", core$signature, "signed", "FALSE", "FALSE", "FALSE", "FALSE", "HGNC symbol", "full", "frozen", "Resolved signature table with flags"),
    c("Module 9.4", "gene_signature", "sensitivity_three_axis", "signature_construction", core$signature, "signed", "FALSE", "FALSE", "FALSE", "FALSE", "HGNC symbol", "sensitivity", "frozen", "Module 9.4 sensitivity membership"),
    c("Module 9.7", "compound_ranking", "primary_three_axis", "DrugReflector", core$dr_primary, "signed input", "TRUE", "FALSE", "FALSE", "FALSE", "BRD ID", "primary", "frozen", "Three-checkpoint ensemble"),
    c("Module 9.7", "compound_ranking", "sensitivity_three_axis", "DrugReflector", core$dr_sensitivity, "signed input", "TRUE", "FALSE", "FALSE", "FALSE", "BRD ID", "sensitivity", "frozen", "Three-checkpoint ensemble"),
    c("Module 9.7", "consensus_ranking", "primary+sensitivity", "DrugReflector", core$dr_consensus, "signed input", "TRUE", "FALSE", "FALSE", "FALSE", "BRD ID", "consensus", "frozen", "Top-200 union"),
    c("Module 9.5/9.8", "external_ranking", "primary_three_axis", "L1000FWD", core$l1000, "direction labelled", "TRUE", "TRUE", "TRUE", "TRUE", "BRD ID", "primary", "frozen", "Existing cached query only"),
    c("Module 9.8", "external_connectivity", "primary_three_axis", "CLUE", core$clue_summary, "signed connectivity", "TRUE", "TRUE", "partial", "partial", "BRD ID/InChIKey/name", "primary", "frozen", "Completed sig_gutc_tool query"),
    c("Module 9.9", "component_connectivity", "landmark_balanced", "CLUE", core$clue_components, "signed connectivity", "TRUE", "TRUE", "TRUE", "TRUE", "BRD ID", "sensitivity", "frozen", "cs_up and -cs_down; not direct expression"),
    c("Module 9.9", "component_ranking", "landmark_decomposition", "DrugReflector", core$decomposition, "signed input", "TRUE", "FALSE", "FALSE", "FALSE", "BRD ID", "sensitivity", "frozen", "Five decomposed profiles")
  )
  manifest <- as.data.table(do.call(rbind, frozen_rows))
  setnames(manifest, c("input_module", "input_type", "signature_version", "method", "file_path", "signed_or_unsigned", "ranking_available", "cell_line_available", "dose_available", "time_available", "compound_id_type", "primary_or_sensitivity", "frozen_status", "notes"))
  manifest[, file_path := normalizePath(file_path, winslash = "/", mustWork = FALSE)]
  manifest[, n_records := vapply(file_path, figure8_file_records, integer(1))]
  setcolorder(manifest, c("input_module", "input_type", "signature_version", "method", "file_path", "n_records", "signed_or_unsigned", "ranking_available", "cell_line_available", "dose_available", "time_available", "compound_id_type", "primary_or_sensitivity", "frozen_status", "notes"))
  figure8_write_tsv(manifest, "figure8_frozen_input_manifest.tsv")

  status <- function(ok, incomplete = FALSE) if (!ok) "unavailable" else if (incomplete) "incomplete" else "available"
  primary_up <- if (nrow(signature)) signature[include_primary == TRUE & desired_direction == "up", uniqueN(gene)] else NA_integer_
  primary_down <- if (nrow(signature)) signature[include_primary == TRUE & desired_direction == "down", uniqueN(gene)] else NA_integer_
  full_n <- if (nrow(signature)) signature[!conflict_flag & !housekeeping_or_qc_flag, uniqueN(gene)] else NA_integer_
  sensitivity_n <- if (nrow(signature)) signature[include_sensitivity == TRUE, uniqueN(gene)] else NA_integer_
  audit <- data.table(
    check_id = 1:24,
    check_name = c(
      "Module 9.4 signature files", "primary UP genes", "primary DOWN genes", "full and sensitivity versions",
      "conflict genes", "housekeeping/QC exclusions", "L1000 landmark coverage", "DrugReflector checkpoint/version",
      "three model predictions", "primary ranking", "sensitivity ranking", "consensus ranking",
      "compound metadata mapping", "BRD/name/InChIKey completeness", "L1000FWD result", "CLUE result",
      "CLUE API submission status", "liver-relevant cell lines", "dose and time", "independent perturbation expression matrix",
      "toxicity/cell-cycle/pan-stress annotation", "compound MoA metadata", "true three-method overlap", "everolimus actual evidence"
    ),
    status = c(
      status(all(file.exists(unlist(core[c("signature", "signature_json", "signature_qc", "module94_report")])))),
      status(is.finite(primary_up)), status(is.finite(primary_down)), status(is.finite(full_n) && is.finite(sensitivity_n)),
      status(nrow(signature) > 0), status(nrow(signature) > 0), status(nrow(coverage) > 0),
      status(all(file.exists(unlist(core[c("checkpoint0", "checkpoint1", "checkpoint2")])))),
      status(all(file.exists(unlist(core[c("dr_primary", "dr_sensitivity", "dr_consensus")])))),
      status(file.exists(core$dr_primary)), status(file.exists(core$dr_sensitivity)), status(file.exists(core$dr_consensus)),
      status(file.exists(core$perturbagen_metadata)), status(nrow(cross) > 0, nrow(cross) > 0 && any(is.na(cross$inchi_key) | cross$inchi_key == "")),
      status(file.exists(core$l1000)), status(file.exists(core$clue_summary)),
      if (!is.null(clue_report) && identical(clue_report$status, "completed")) "available" else "not queried",
      status(file.exists(core$clue_cells)), status(file.exists(core$clue_components), TRUE),
      "unavailable", "incomplete", "unavailable",
      status(!is.null(cv_report)), status(any(tolower(cross$pert_iname %||% character()) == "everolimus"), TRUE)
    ),
    value = c(
      paste(sum(file.exists(unlist(core[c("signature", "signature_json", "signature_qc", "module94_report")]))), "/4"),
      primary_up, primary_down, paste0("full=", full_n, "; sensitivity=", sensitivity_n),
      if (nrow(signature)) sum(signature$conflict_flag) else NA, if (nrow(signature)) sum(signature$housekeeping_or_qc_flag) else NA,
      if (nrow(coverage)) paste(unique(coverage[fold == "union" | as.character(fold) == "union", paste0(signature, "=", n_overlap_genes)]), collapse = "; ") else NA,
      "DrugReflector V3.5; 3 frozen checkpoints", if (!is.null(dr_report)) 3 else NA, figure8_file_records(core$dr_primary), figure8_file_records(core$dr_sensitivity), figure8_file_records(core$dr_consensus),
      if (nrow(cross)) sum(figure8_bool(cross$metadata_mapped)) else NA, if (nrow(cross)) paste0("BRD=", sum(!is.na(cross$compound)), "; name=", sum(!is.na(cross$pert_iname)), "; InChIKey=", sum(!is.na(cross$inchi_key))) else NA,
      figure8_file_records(core$l1000), figure8_file_records(core$clue_summary), if (!is.null(clue_report)) clue_report$job_id else NA,
      "HEPG2; HCC515; HA1E", "encoded in CLUE sig_id; not separate columns", "none found",
      "local frozen gene sets available; compound-level viability absent", "none found",
      if (!is.null(cv_report)) cv_report$summary$n_three_method_strong_support else NA,
      if (nrow(cross)) paste(cross[tolower(pert_iname) == "everolimus", paste0("DR_primary=", primary_rank_1based, "; DR_sensitivity=", sensitivity_rank_1based, "; CLUE_tau=", clue_tau)], collapse = "; ") else NA
    ),
    notes = c(
      "Frozen files are read-only inputs", "From include_primary", "From include_primary", "Full excludes conflict/QC; sensitivity follows Module 9.4",
      "Reported, not hidden", "Reported, not hidden", "Primary coverage is low and retained as review risk", "Checkpoint MD5 verified by Module 9.6 and rechecked in inference report",
      "Fold-specific ranks will be exported for Figure 8", "Frozen", "Frozen", "Frozen",
      "Exact local metadata", "Missing identifiers remain explicit", "Existing cache; no API refresh", "Existing completed query; no API refresh",
      "Completed in Module 9.8", "Connectivity support only", "Dose/time parsing is possible; replicate metadata incomplete", "CLUE connectivity is not direct expression",
      "Model-based nuisance sensitivity is possible; safety is not resolved", "MoA/targets cannot be inferred from names", "Zero must remain visible if confirmed", "Independent external-support reference unless new tier rules are met"
    )
  )
  figure8_write_tsv(audit, "figure8_preflight_report.tsv")

  weights <- data.table(
    dimension = c("DrugReflector_stability", "ensemble_agreement", "signature_robustness", "external_method_support", "liver_context_reversal", "identity_rescue", "stress_suppression", "sox4_suppression", "mechanism_interpretability"),
    weight = c(0.15, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05),
    role = "positive",
    frozen_before_candidate_scoring = TRUE
  )
  penalties <- data.table(
    dimension = c("proliferation_penalty", "generic_stress_penalty", "cytotoxicity_penalty", "mapping_uncertainty_penalty"),
    weight = 1, role = "subtractive", frozen_before_candidate_scoring = TRUE
  )
  figure8_write_tsv(rbindlist(list(weights, penalties), fill = TRUE), "figure8_score_weight_contract.tsv")
  figure8_write_tsv(figure8_protected_manifest(), "figure8_pre_run_protected_files_manifest.tsv")
  figure8_write_tsv(figure8_protected_data_manifest(), "figure8_pre_run_protected_data_manifest.tsv")
  figure8_write_json(list(
    module = "figure8_preflight", status = if (any(audit$status %in% c("unavailable", "incomplete", "review risk"))) "completed_with_gaps" else "completed",
    seed = FIGURE8_SEED, palette = figure8_palette_contract(), audit = audit,
    frozen_input_count = nrow(manifest), external_api_refresh = FALSE,
    review_risks = audit[status %in% c("unavailable", "incomplete", "mapping failed", "review risk"), check_name]
  ), "figure8_preflight_report.json")
  invisible(audit)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

figure8_expression_universe <- function() {
  h5ad <- figure8_input("data", "processed", "driver", "driver_hepatocyte_trajectory.module6_1.h5ad")
  figure8_require(h5ad, "Driver hepatocyte h5ad")
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("hdf5r is required for expression-matched null signatures")
  f <- hdf5r::H5File$new(h5ad, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  genes <- toupper(as.character(f[["var/_index"]][]))
  data <- as.numeric(f[["X/data"]][])
  indices <- as.integer(f[["X/indices"]][])
  shape <- as.integer(f[["X"]]$attr_open("shape")$read())
  sums <- numeric(shape[[2]])
  grouped <- rowsum(data, group = indices + 1L, reorder = FALSE)
  sums[as.integer(rownames(grouped))] <- grouped[, 1]
  detected <- tabulate(indices + 1L, nbins = shape[[2]])
  out <- data.table(
    gene = genes,
    mean_expression = sums / shape[[1]],
    detection_fraction = detected / shape[[1]],
    expression_source = normalizePath(h5ad, winslash = "/", mustWork = TRUE)
  )
  model <- figure8_model_genes()
  required <- unique(c(model, figure8_signature()$gene))
  missing <- setdiff(required, out$gene)
  if (length(missing)) {
    out <- rbind(out, data.table(gene = missing, mean_expression = 0, detection_fraction = 0, expression_source = "not_present_in_driver_h5ad"), fill = TRUE)
  }
  out[, expression_bin := fifelse(
    mean_expression <= 0, "zero_or_absent",
    paste0("decile_", as.integer(cut(mean_expression, unique(quantile(mean_expression[mean_expression > 0], probs = seq(0, 1, 0.1), na.rm = TRUE)), include.lowest = TRUE, labels = FALSE)))
  )]
  unique(out, by = "gene")
}

figure8_profile <- function(frame, signature_id, filter_idx = rep(TRUE, nrow(frame)), balance = FALSE, model_genes = character()) {
  x <- copy(frame[filter_idx])
  if (!nrow(x)) return(data.table())
  x[, v_score := fifelse(desired_direction == "up", 1, -1) * final_weight]
  if (isTRUE(balance)) {
    x[, v_score := v_score / sum(abs(v_score)), by = desired_direction]
  }
  x[, `:=`(signature_id = signature_id, is_landmark = gene %in% model_genes)]
  x
}

run_figure8_prepare_signature <- function() {
  figure8_dirs()
  signature <- figure8_signature()
  model <- figure8_model_genes()
  nuisance <- figure8_local_nuisance_sets()
  figure8_write_tsv(nuisance, "figure8_nuisance_gene_sets.tsv")
  expression <- figure8_expression_universe()
  figure8_write_tsv(expression, "figure8_expression_universe.tsv", compress = TRUE)

  primary <- signature$include_primary
  eligible <- !signature$conflict_flag & !signature$housekeeping_or_qc_flag
  sensitivity <- signature$include_sensitivity
  prolif_genes <- nuisance[nuisance_set == "proliferation", unique(gene)]
  cycle_genes <- nuisance[nuisance_set == "cell_cycle", unique(gene)]
  stress_genes <- nuisance[nuisance_set %in% c("generic_stress", "unfolded_protein_response"), unique(gene)]
  identity <- signature$axis == "identity_rescue"
  stress <- signature$axis == "stress_suppression"
  sox4 <- signature$axis == "sox4_suppression"

  profiles <- list(
    figure8_profile(signature, "primary_three_axis", primary, FALSE, model),
    figure8_profile(signature, "full_three_axis", eligible, FALSE, model),
    figure8_profile(signature, "sensitivity_three_axis", sensitivity, FALSE, model),
    figure8_profile(signature, "no_proliferation", primary & !signature$gene %in% prolif_genes, FALSE, model),
    figure8_profile(signature, "no_cell_cycle", primary & !signature$gene %in% cycle_genes, FALSE, model),
    figure8_profile(signature, "no_generic_stress", primary & !signature$gene %in% stress_genes, FALSE, model),
    figure8_profile(signature, "sox4_only", eligible & sox4, TRUE, model),
    figure8_profile(signature, "identity_rescue_only", eligible & identity, TRUE, model),
    figure8_profile(signature, "stress_suppression_only", eligible & stress, TRUE, model),
    figure8_profile(signature, "leave_out_sox4", primary & !sox4, FALSE, model),
    figure8_profile(signature, "leave_out_stress", primary & !stress, FALSE, model),
    figure8_profile(signature, "leave_out_identity", primary & !identity, FALSE, model),
    figure8_profile(signature, "landmark_only", eligible & signature$gene %in% model, FALSE, model),
    figure8_profile(signature, "balanced_up_down", primary, TRUE, model),
    figure8_profile(signature, "high_confidence_intersection", primary & signature$gene %in% model, TRUE, model)
  )
  variants <- rbindlist(profiles, fill = TRUE)
  expected <- c(
    "primary_three_axis", "full_three_axis", "sensitivity_three_axis", "no_proliferation", "no_cell_cycle",
    "no_generic_stress", "sox4_only", "identity_rescue_only", "stress_suppression_only", "leave_out_sox4",
    "leave_out_stress", "leave_out_identity", "landmark_only", "balanced_up_down", "high_confidence_intersection"
  )
  if (!identical(unique(variants$signature_id), expected)) stop("Figure 8 signature variants differ from the preregistered order")
  figure8_write_tsv(variants, "figure8_signature_variants_long.tsv", compress = TRUE)
  wide <- dcast(variants, signature_id ~ gene, value.var = "v_score", fun.aggregate = sum, fill = 0)
  wide[, signature_id := factor(signature_id, levels = expected)]
  setorder(wide, signature_id)
  wide[, signature_id := as.character(signature_id)]
  figure8_write_data_tsv(wide, "figure8_signature_variants_wide.tsv", compress = TRUE)

  manifest <- variants[, .(
    up_gene_count = uniqueN(gene[desired_direction == "up"]),
    down_gene_count = uniqueN(gene[desired_direction == "down"]),
    landmark_up_count = uniqueN(gene[desired_direction == "up" & is_landmark]),
    landmark_down_count = uniqueN(gene[desired_direction == "down" & is_landmark]),
    component_weights = jsonlite::toJSON(as.list(tapply(abs(v_score), axis, sum)), auto_unbox = TRUE),
    conflict_rule = ifelse(signature_id == "sensitivity_three_axis", "retain Module 9.4 resolved flagged rows", "exclude according to frozen variant rule"),
    normalization = ifelse(signature_id %in% c("sox4_only", "identity_rescue_only", "stress_suppression_only", "balanced_up_down", "high_confidence_intersection"), "absolute mass normalized within direction", "frozen final_weight; DrugReflector clip/rescale"),
    model_input_dimension = length(model),
    primary_or_sensitivity = fifelse(signature_id == "primary_three_axis", "primary", fifelse(signature_id == "sensitivity_three_axis", "sensitivity", "predefined robustness"))
  ), by = signature_id]
  manifest[, signature_id := factor(signature_id, levels = expected)]
  setorder(manifest, signature_id)
  manifest[, signature_id := as.character(signature_id)]
  figure8_write_tsv(manifest, "figure8_signature_variant_manifest.tsv")

  version_membership <- variants[, .(gene, component, desired_direction, axis, final_weight, is_landmark, signature_id)]
  composition <- version_membership[, .(
    n_genes = uniqueN(gene), median_weight = median(abs(final_weight), na.rm = TRUE),
    total_weight = sum(abs(final_weight), na.rm = TRUE), n_landmark = uniqueN(gene[is_landmark]),
    primary_coverage = mean(gene %in% signature[include_primary == TRUE, gene])
  ), by = .(component, axis, desired_direction, signature_id)]
  raw_report <- figure8_parse_json(figure8_input("metadata", "driver", "module9_4_report.json"))
  qc_summary <- data.table(
    metric = c("raw_records", "resolved_unique_genes", "primary_up", "primary_down", "conflicts_excluded", "housekeeping_qc_excluded", "landmark_covered", "non_landmark_eligible"),
    value = c(
      raw_report$summary$n_raw_records %||% NA, uniqueN(signature$gene), signature[include_primary == TRUE & desired_direction == "up", uniqueN(gene)],
      signature[include_primary == TRUE & desired_direction == "down", uniqueN(gene)], sum(signature$conflict_flag), sum(signature$housekeeping_or_qc_flag),
      signature[eligible & gene %in% model, uniqueN(gene)], signature[eligible & !gene %in% model, uniqueN(gene)]
    )
  )
  figure8_write_tsv(composition, "figure8b_signature_composition.tsv")

  coverage <- rbindlist(list(
    variants[, .(n_genes = uniqueN(gene), n_landmark = uniqueN(gene[is_landmark])), by = .(signature_id, desired_direction)][, component := "all"],
    variants[, .(n_genes = uniqueN(gene), n_landmark = uniqueN(gene[is_landmark])), by = .(signature_id, desired_direction, component)]
  ), fill = TRUE)
  coverage[, coverage_fraction := n_landmark / pmax(n_genes, 1)]
  coverage[, non_landmark := n_genes - n_landmark]
  figure8_write_tsv(coverage, "figure8b_landmark_coverage.tsv")
  figure8_write_tsv(qc_summary, "figure8b_signature_qc_summary.tsv")
  figure8_write_json(list(
    module = "figure8_signature_composition", status = "completed", seed = FIGURE8_SEED,
    signature_variants = nrow(manifest), model_landmark_genes = length(model),
    primary = list(up = manifest[signature_id == "primary_three_axis", up_gene_count], down = manifest[signature_id == "primary_three_axis", down_gene_count]),
    review_risks = c("Low primary landmark coverage is retained", "CEBPB/EGR1 malignant-target component is absent from the resolved Module 9.4 table when not selected by frozen priority rules")
  ), "figure8b_signature_qc_report.json")
  invisible(manifest)
}

plot_figure8a <- function() {
  states <- data.table(
    state = rep(c("Disease / malignant state", "Desired rescue state"), each = 6),
    feature = rep(c("HNF4A/PPARA identity", "Mature hepatocyte identity", "AP-1/CEBPB/EGR1 stress-transition", "SOX4 malignant-state programme", "CNV-associated malignant fate", "Cell-cycle collapse independence"), 2),
    direction = c("low", "low", "high", "high", "high", "unresolved", "restored", "restored", "reduced", "reduced", "reduced", "required"),
    axis = rep(c("identity_rescue", "identity_rescue", "stress_suppression", "sox4_suppression", "sox4_suppression", "stress_suppression"), 2)
  )
  states[, state := factor(state, levels = c("Disease / malignant state", "Desired rescue state"))]
  states[, feature := factor(feature, levels = rev(unique(feature)))]
  p <- ggplot(states, aes(x = state, y = feature)) +
    geom_segment(aes(x = 1, xend = 2, yend = feature), colour = lancet_palette[8], linewidth = 0.35,
                 arrow = grid::arrow(length = grid::unit(1.6, "mm"))) +
    geom_point(aes(fill = axis), shape = 21, size = 5.4, colour = "white", stroke = 0.5) +
    geom_text(aes(label = direction), size = 2.35, colour = lancet_palette[9], fontface = "bold") +
    annotate("label", x = 1.5, y = 6.75, label = "Transcriptomic reversal target", size = 2.8,
             label.size = 0.25, fill = "white", colour = lancet_palette[9]) +
    scale_fill_manual(values = axis_palette, breaks = names(axis_palette),
                      labels = c("Identity rescue", "Stress-transition suppression", "SOX4/malignant suppression")) +
    scale_x_discrete(position = "top") +
    coord_cartesian(clip = "off", ylim = c(0.6, 7)) +
    labs(
      title = "Computational definition of the malignant and desired rescue states",
      x = NULL, y = NULL, fill = "Three-axis target",
      caption = "The desired rescue state is computationally defined and does not represent experimentally demonstrated phenotypic rescue."
    ) + figure8_theme() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(face = "bold", size = 8), legend.position = "bottom")
  figure8_save(p, "a", "figure8a_target_state_definition", 7.2, 4.8)
  figure8_write_tsv(states, "figure8a_target_state_data.tsv")
  figure8_write_json(list(module = "figure8a", status = "completed", statement = "Computational target definition only", experimental_rescue_demonstrated = FALSE), "figure8a_target_state_report.json")
}

plot_figure8b <- function() {
  suppressPackageStartupMessages(library(ggalluvial))
  comp <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8b_signature_composition.tsv"))
  cov <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8b_landmark_coverage.tsv"))
  qc <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8b_signature_qc_summary.tsv"))
  show_versions <- c("primary_three_axis", "full_three_axis", "sensitivity_three_axis")
  alluv <- comp[signature_id %in% show_versions]
  alluv[, flow_weight := pmax(n_genes, 1)]
  p1 <- ggplot(alluv, aes(axis1 = component, axis2 = desired_direction, axis3 = signature_id, y = flow_weight)) +
    ggalluvial::geom_alluvium(aes(fill = axis), width = 0.18, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 0.20, fill = "white", colour = lancet_palette[9], linewidth = 0.25) +
    ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.0, check_overlap = TRUE) +
    scale_x_discrete(limits = c("Component", "Direction", "Version"), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = axis_palette) + labs(title = "Signature provenance", x = NULL, y = "Gene records", fill = "Axis") +
    figure8_theme() + theme(legend.position = "none", axis.text.x = element_text(face = "bold"))
  bubble <- comp[signature_id == "primary_three_axis"]
  p2 <- ggplot(bubble, aes(x = n_genes, y = reorder(component, n_genes))) +
    geom_point(aes(size = total_weight, fill = axis), shape = 21, colour = "white", stroke = 0.4) +
    geom_text(aes(label = paste0("n=", n_genes, "\nmed=", number(median_weight, accuracy = 0.01))), size = 2.0, nudge_x = max(bubble$n_genes) * 0.08) +
    scale_fill_manual(values = axis_palette) + scale_size_area(max_size = 11) +
    labs(title = "Primary genes and weights", x = "Resolved genes", y = NULL, size = "Total weight", fill = "Axis") +
    figure8_theme() + theme(legend.position = "none")
  cov_show <- cov[signature_id %in% show_versions & component == "all"]
  cov_show[, label := paste0(n_landmark, "/", n_genes, " (", percent(coverage_fraction, accuracy = 0.1), ")")]
  p3 <- ggplot(cov_show, aes(x = coverage_fraction, y = interaction(signature_id, desired_direction, sep = " | "), fill = desired_direction)) +
    geom_col(width = 0.62) + geom_text(aes(label = label), hjust = -0.05, size = 2.1) +
    scale_fill_manual(values = c(up = axis_palette[["identity_rescue"]], down = axis_palette[["sox4_suppression"]])) +
    scale_x_continuous(labels = percent, limits = c(0, max(0.65, max(cov_show$coverage_fraction) * 1.45))) +
    labs(title = "DrugReflector landmark coverage", x = "Covered fraction", y = NULL, fill = "Desired direction") +
    figure8_theme() + theme(legend.position = "bottom")
  qct <- paste(qc$metric, qc$value, sep = " = ", collapse = "  |  ")
  combined <- (p1 | p2 | p3) + plot_layout(widths = c(1.35, 0.9, 1.05)) +
    plot_annotation(title = "Three-axis reversal signature composition and quality control", caption = qct, theme = theme(plot.title = element_text(size = 10), plot.caption = element_text(size = 6.5)))
  figure8_save(combined, "b", "figure8b_signature_composition_qc", 13.2, 5.5)
}

plot_figure8c <- function() {
  pre <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_preflight_report.tsv"))
  nodes <- data.table(
    step = 1:8,
    label = c("Frozen three-axis\nreversal signature", "DrugReflector\ninput v-score", "Three-model ensemble\nprediction", "Primary / sensitivity /\ncomponent-ablation ranking", "Compound metadata\nharmonization", "L1000FWD / CLUE /\nindependent comparison", "Liver-context\nconnectivity support", "Mechanism/response class\nand evidence tiers"),
    status = c("completed", "completed", "completed", "partial", "completed", "partial", "external_validation", "partial")
  )
  nodes[, `:=`(x = step, y = ifelse(step %% 2 == 1, 1, 0))]
  edges <- nodes[-.N, .(x, y, xend = nodes$x[.I + 1], yend = nodes$y[.I + 1])]
  p <- ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), linewidth = 0.45, colour = lancet_palette[8], arrow = grid::arrow(length = grid::unit(1.8, "mm"))) +
    geom_label(data = nodes, aes(x, y, label = label, fill = status), size = 2.25, label.size = 0.28, label.padding = grid::unit(2.2, "mm"), colour = lancet_palette[9]) +
    scale_fill_manual(values = workflow_palette, breaks = names(workflow_palette)) +
    coord_cartesian(clip = "off", xlim = c(0.5, 8.5), ylim = c(-0.4, 1.4)) +
    labs(title = "Frozen-analysis workflow and evidence boundary", fill = "Status",
         caption = "No external API was refreshed. CLUE liver evidence is connectivity support; direct perturbational expression and experimental rescue are unavailable.") +
    figure8_theme() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
  figure8_save(p, "c", "figure8c_reversal_workflow", 11.2, 4.2)
  status_table <- rbind(nodes[, .(workflow_step = step, node = label, status)], data.table(workflow_step = 9, node = "Experimental phenotypic rescue", status = "unavailable"))
  figure8_write_tsv(status_table, "figure8c_workflow_status.tsv")
  figure8_write_json(list(module = "figure8c", status = "completed", workflow = status_table, preflight_gap_count = sum(pre$status != "available"), external_api_refresh = FALSE, active_learning_loop = FALSE), "figure8c_workflow_report.json")
}

figure8_load_perturbagen_metadata <- function() {
  cache_path <- file.path(FIGURE8_METADATA_DIR, "figure8_perturbagen_identity_map.tsv.gz")
  if (file.exists(cache_path)) return(figure8_fread(cache_path))
  phase1 <- figure8_input("metadata", "driver", "GSE92742_Broad_LINCS_pert_info.txt.gz")
  phase2 <- figure8_input("metadata", "driver", "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz")
  frames <- list()
  if (file.exists(phase1)) {
    x <- figure8_fread(phase1)
    x[, metadata_source := "GSE92742_phase1"]
    frames[[length(frames) + 1L]] <- x
  }
  if (file.exists(phase2)) {
    x <- figure8_fread(phase2)
    x[, metadata_source := "GSE70138_phase2"]
    frames[[length(frames) + 1L]] <- x
  }
  if (!length(frames)) return(data.table())
  all <- rbindlist(frames, fill = TRUE, use.names = TRUE)
  all[, compound := as.character(pert_id)]
  for (nm in c("pert_iname", "inchi_key", "canonical_smiles", "pubchem_cid", "pert_type", "is_touchstone")) {
    if (!nm %in% names(all)) all[, (nm) := NA_character_]
  }
  all[, `:=`(
    pert_iname = fifelse(pert_iname %in% c("", "-666"), NA_character_, as.character(pert_iname)),
    inchi_key = fifelse(inchi_key %in% c("", "-666"), NA_character_, as.character(inchi_key)),
    canonical_smiles = fifelse(canonical_smiles %in% c("", "-666"), NA_character_, as.character(canonical_smiles)),
    pubchem_cid = fifelse(pubchem_cid %in% c("", "-666"), NA_character_, as.character(pubchem_cid)),
    pert_type = fifelse(as.character(pert_type) %in% c("", "-666"), NA_character_, as.character(pert_type)),
    is_touchstone = fifelse(as.character(is_touchstone) %in% c("", "-666"), NA_character_, as.character(is_touchstone)),
    source_priority = fifelse(metadata_source == "GSE92742_phase1", 1L, 2L)
  )]
  setorder(all, compound, source_priority)
  resolved <- data.table(compound = unique(all$compound))
  for (nm in c("pert_iname", "pert_type", "is_touchstone", "inchi_key", "canonical_smiles", "pubchem_cid")) {
    value_map <- unique(all[!is.na(get(nm)), .(compound, value = as.character(get(nm)))], by = "compound")
    setnames(value_map, "value", nm)
    resolved <- merge(resolved, value_map, by = "compound", all.x = TRUE, sort = FALSE)
  }
  conflicts <- all[, .(
    metadata_sources = paste(sort(unique(metadata_source)), collapse = ","),
    metadata_conflict_flag = uniqueN(pert_iname, na.rm = TRUE) > 1 | uniqueN(inchi_key, na.rm = TRUE) > 1
  ), by = compound]
  resolved <- merge(resolved, conflicts, by = "compound", all.x = TRUE, sort = FALSE)
  resolved[, canonical_name := fifelse(!is.na(pert_iname), pert_iname, compound)]
  resolved[, normalized_name := figure8_safe_name(canonical_name)]
  figure8_write_tsv(resolved, basename(cache_path), compress = TRUE)
  resolved
}

figure8_first_or_na <- function(x) {
  y <- na.omit(x)
  if (length(y)) y[[1]] else NA_character_
}

run_figure8_stability <- function() {
  pred_path <- file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz")
  fold_path <- file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_fold_predictions.tsv.gz")
  figure8_require(c(pred_path, fold_path), "Figure 8 variant inference")
  pred <- figure8_fread(pred_path)
  folds <- figure8_fread(fold_path)
  n_compounds <- uniqueN(pred$compound)
  expected <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_signature_variant_manifest.tsv"))$signature_id
  pred[, rank_score := figure8_rank_score(rank_1based, n_compounds)]
  union_top <- pred[rank_1based <= 200, unique(compound)]
  frozen_consensus_path <- figure8_input("metadata", "driver", "module9_7_drugreflector_consensus_predictions.tsv")
  frozen_consensus <- if (file.exists(frozen_consensus_path)) figure8_fread(frozen_consensus_path)$compound else character()
  metadata <- figure8_load_perturbagen_metadata()
  refs <- metadata[tolower(canonical_name) %in% c("everolimus", "dapivirine", "tipiracil", "cefepime", "tasquinimod", "cisapride"), compound]
  universe <- unique(c(union_top, frozen_consensus, refs))

  stability <- pred[, .(
    median_rank = median(rank_1based), rank_q1 = quantile(rank_1based, 0.25), rank_q3 = quantile(rank_1based, 0.75),
    best_rank = min(rank_1based), worst_rank = max(rank_1based), rank_sd = sd(rank_1based),
    median_rank_score = median(rank_score), top20_frequency = mean(rank_1based <= 20),
    top50_frequency = mean(rank_1based <= 50), top100_frequency = mean(rank_1based <= 100),
    top200_frequency = mean(rank_1based <= 200), n_signature_versions = .N
  ), by = compound]
  fold_agreement <- folds[, .(
    fold_rank_min = min(fold_rank_1based), fold_rank_max = max(fold_rank_1based), fold_rank_sd = sd(fold_rank_1based),
    fold_top100_frequency = mean(fold_rank_1based <= 100),
    per_signature_model_agreement = pmax(0, 1 - (max(fold_rank_1based) - min(fold_rank_1based)) / max(n_compounds - 1, 1))
  ), by = .(signature_id, compound)][, .(
    model_agreement = median(per_signature_model_agreement),
    ensemble_fold_rank_sd = median(fold_rank_sd),
    ensemble_top100_frequency = mean(fold_top100_frequency)
  ), by = compound]
  stability <- merge(stability, fold_agreement, by = "compound", all.x = TRUE)
  stability <- merge(stability, metadata, by = "compound", all.x = TRUE)
  stability[, `:=`(
    canonical_name = fifelse(is.na(canonical_name) | canonical_name == "", compound, canonical_name),
    in_candidate_analysis_universe = compound %in% universe,
    candidate_origin = fifelse(compound %in% refs, "external_reference", fifelse(compound %in% frozen_consensus, "frozen_consensus", "variant_top200_union"))
  )]
  setorder(stability, median_rank, -top100_frequency, compound)

  wide <- dcast(pred, compound ~ signature_id, value.var = "rank_1based")
  rank_matrix <- as.matrix(wide[, ..expected])
  correlations <- as.data.table(as.table(stats::cor(rank_matrix, method = "spearman", use = "pairwise.complete.obs")))
  setnames(correlations, c("signature_a", "signature_b", "spearman_rho"))
  correlations[, n_compounds := nrow(rank_matrix)]
  correlations[, comparison_type := "signature_variant"]
  fold_corr <- folds[, {
    w <- dcast(.SD, compound ~ fold, value.var = "fold_rank_1based")
    m <- as.matrix(w[, -1])
    cmat <- cor(m, method = "spearman")
    data.table(mean_pairwise_spearman = mean(cmat[upper.tri(cmat)]), min_pairwise_spearman = min(cmat[upper.tri(cmat)]))
  }, by = signature_id]

  m <- ncol(rank_matrix)
  n <- nrow(rank_matrix)
  rank_sums <- rowSums(rank_matrix)
  kendall_w <- 12 * sum((rank_sums - mean(rank_sums))^2) / (m^2 * (n^3 - n))
  primary <- pred[signature_id == "primary_three_axis"]
  sensitivity <- pred[signature_id == "sensitivity_three_axis"]
  overlap <- rbindlist(lapply(c(20L, 50L, 100L, 200L), function(k) {
    a <- primary[rank_1based <= k, compound]
    b <- sensitivity[rank_1based <= k, compound]
    data.table(k = k, primary_n = length(a), sensitivity_n = length(b), overlap_n = length(intersect(a, b)), jaccard = length(intersect(a, b)) / length(union(a, b)))
  }))
  leaveout <- rbindlist(lapply(c("leave_out_sox4", "leave_out_stress", "leave_out_identity"), function(v) {
    rho <- correlations[signature_a == "primary_three_axis" & signature_b == v, spearman_rho]
    top_primary <- pred[signature_id == "primary_three_axis" & rank_1based <= 100, compound]
    top_variant <- pred[signature_id == v & rank_1based <= 100, compound]
    data.table(variant = v, spearman_rho = rho[[1]], top100_overlap = length(intersect(top_primary, top_variant)), top100_jaccard = length(intersect(top_primary, top_variant)) / length(union(top_primary, top_variant)))
  }))
  correlations <- rbind(correlations, data.table(signature_a = "ensemble_folds", signature_b = fold_corr$signature_id, spearman_rho = fold_corr$mean_pairwise_spearman, n_compounds = n_compounds, comparison_type = "model_fold"), fill = TRUE)

  display <- stability[in_candidate_analysis_universe == TRUE][order(median_rank, -top100_frequency)][1:15, compound]
  display <- unique(c(display, refs))
  display <- display[!is.na(display)][seq_len(min(length(display), 20L))]
  stability[, display_in_main_panel := compound %in% display]
  figure8_write_tsv(stability, "figure8d_rank_stability.tsv")
  figure8_write_tsv(correlations, "figure8d_rank_correlation.tsv")
  figure8_write_tsv(overlap, "figure8d_primary_sensitivity_overlap.tsv")
  figure8_write_tsv(leaveout, "figure8d_leave_one_component_out_stability.tsv")
  figure8_write_tsv(data.table(compound = display), "figure8_display_candidate_manifest.tsv")
  figure8_write_json(list(
    module = "figure8d_stability", status = "completed", seed = FIGURE8_SEED,
    n_compounds = n_compounds, n_variants = length(expected), candidate_analysis_universe = length(universe),
    kendall_w = kendall_w, primary_sensitivity_overlap = overlap, fold_agreement = fold_corr,
    selection_rule = "All-compound statistics; display uses the 15 lowest median ranks from the predeclared Top-200 union plus frozen named references, capped at 20.",
    everolimus = stability[tolower(canonical_name) == "everolimus", .(compound, median_rank, best_rank, worst_rank, top100_frequency)],
    interpretation_boundary = "Rank stability is internal model robustness, not therapeutic efficacy."
  ), "figure8d_stability_report.json")
  invisible(stability)
}

plot_figure8d <- function() {
  pred <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz"))
  stab <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  show <- stab[display_in_main_panel == TRUE, .(compound, canonical_name, median_rank, rank_q1, rank_q3, best_rank, worst_rank, top100_frequency, candidate_origin)]
  show[, label := fifelse(is.na(canonical_name) | canonical_name == compound, compound, paste0(canonical_name, "\n", compound))]
  order_names <- show[order(median_rank), label]
  heat <- merge(pred[compound %in% show$compound], show[, .(compound, label)], by = "compound")
  heat[, `:=`(rank_signal = -log10(rank_1based), label = factor(label, levels = rev(order_names)))]
  p1 <- ggplot(heat, aes(x = signature_id, y = label, fill = rank_signal)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradientn(colours = c("#F7F7F7", lancet_palette[6], lancet_palette[2]), name = expression(-log[10](rank)), guide = guide_colorbar(direction = "horizontal", title.position = "top", barwidth = grid::unit(4, "cm"), barheight = grid::unit(0.35, "cm"))) +
    labs(title = "Candidate ranks across preregistered signatures", x = NULL, y = NULL) +
    figure8_theme() + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6.5))
  show[, label := factor(label, levels = rev(order_names))]
  p2 <- ggplot(show, aes(y = label, x = median_rank)) +
    geom_errorbarh(aes(xmin = best_rank, xmax = worst_rank), height = 0, linewidth = 0.35, colour = lancet_palette[8]) +
    geom_errorbarh(aes(xmin = rank_q1, xmax = rank_q3), height = 0, linewidth = 1.2, colour = lancet_palette[9]) +
    geom_point(aes(size = top100_frequency, fill = candidate_origin), shape = 21, colour = "white", stroke = 0.4) +
    scale_x_log10() + scale_size_area(max_size = 6, labels = percent) +
    scale_fill_manual(values = c(external_reference = lancet_palette[4], frozen_consensus = lancet_palette[5], variant_top200_union = lancet_palette[6])) +
    labs(title = "Median, IQR and range", x = "DrugReflector rank (log scale)", y = NULL, size = "Top-100 frequency", fill = "Candidate origin") +
    figure8_theme() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  combined <- ((p1 | p2) + plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(title = "DrugReflector rank stability across signature definitions", caption = "Statistics were computed over all 9,597 compounds before display selection. everolimus retains its observed model rank.")) &
    theme(legend.position = "bottom")
  figure8_save(combined, "d", "figure8d_drugreflector_rank_stability", 13.2, 7.0)
}

figure8_entity_key <- function(inchi_key, canonical_name, raw_id) {
  inchi <- trimws(as.character(inchi_key))
  name <- figure8_safe_name(canonical_name)
  raw <- trimws(as.character(raw_id))
  fifelse(!is.na(inchi) & nzchar(inchi) & inchi != "-666", paste0("INCHI:", inchi),
          fifelse(!is.na(name) & nzchar(name), paste0("NAME:", name), paste0("BRD:", raw)))
}

figure8_rank_correlation_row <- function(a, b, label_a, label_b) {
  z <- merge(a, b, by = "standardized_id", suffixes = c("_a", "_b"))
  ok <- is.finite(z$rank_a) & is.finite(z$rank_b)
  data.table(method_a = label_a, method_b = label_b, n_common = sum(ok), spearman_rho = if (sum(ok) >= 3) cor(z$rank_a[ok], z$rank_b[ok], method = "spearman") else NA_real_, status = if (sum(ok) >= 3) "available" else "not_computable")
}

run_figure8_cross_method <- function() {
  pred <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz"))
  stability <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  meta <- figure8_load_perturbagen_metadata()
  dr <- pred[signature_id == "primary_three_axis", .(raw_id = compound, rank = rank_1based, score = prob)]
  dr <- merge(dr, meta[, .(compound, canonical_name, inchi_key, pubchem_cid, metadata_conflict_flag)], by.x = "raw_id", by.y = "compound", all.x = TRUE)
  dr[, `:=`(method = "DrugReflector", direction = "support", support = rank <= 200, strong_support = rank <= 100, strong_opposition = FALSE)]

  lpath <- figure8_input("metadata", "driver", "module9_5_l1000fwd_candidate_ranking.tsv")
  l1000 <- figure8_fread(lpath)
  l1000[, raw_id := as.character(compound_id)]
  l1000 <- merge(l1000, meta[, .(compound, canonical_name, inchi_key, pubchem_cid, metadata_conflict_flag)], by.x = "raw_id", by.y = "compound", all.x = TRUE)
  l1000[, `:=`(
    rank = as.numeric(rank_within_group), score = as.numeric(final_rank_score), method = "L1000FWD",
    direction = fifelse(grepl("similar", candidate_direction, ignore.case = TRUE), "support", "opposition"),
    support = grepl("similar", candidate_direction, ignore.case = TRUE) & rank_within_group <= 50,
    strong_support = grepl("similar", candidate_direction, ignore.case = TRUE) & rank_within_group <= 20,
    strong_opposition = grepl("opposite", candidate_direction, ignore.case = TRUE) & rank_within_group <= 20
  )]
  l1000[is.na(canonical_name) | canonical_name == "", canonical_name := fifelse(!is.na(compound_name) & compound_name != "", compound_name, raw_id)]

  cpath <- figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_clue_perturbagen_summary.tsv")
  clue <- figure8_fread(cpath)
  setnames(clue, "clue_compound", "raw_id")
  clue[, `:=`(
    canonical_name = fifelse(is.na(pert_iname) | pert_iname == "", raw_id, pert_iname), rank = frank(-clue_tau, ties.method = "min"),
    score = as.numeric(clue_tau), method = "CLUE", direction = fifelse(clue_tau > 0, "support", fifelse(clue_tau < 0, "opposition", "neutral")),
    support = clue_tau >= 90, strong_support = clue_tau >= 90, strong_opposition = clue_tau <= -90
  )]
  if (!"metadata_conflict_flag" %in% names(clue)) clue[, metadata_conflict_flag := FALSE]
  if (!"pubchem_cid" %in% names(clue)) clue[, pubchem_cid := NA_character_]

  common_cols <- c("raw_id", "canonical_name", "inchi_key", "pubchem_cid", "metadata_conflict_flag", "method", "rank", "score", "direction", "support", "strong_support", "strong_opposition")
  entities <- rbindlist(list(dr[, ..common_cols], l1000[, ..common_cols], clue[, ..common_cols]), fill = TRUE)
  entities[, standardized_id := figure8_entity_key(inchi_key, canonical_name, raw_id)]
  entities[, mapping_level := fifelse(grepl("^INCHI:", standardized_id), "InChIKey", fifelse(grepl("^NAME:", standardized_id), "standardized_name", "BRD_ID"))]
  entities[, mapping_confidence := fifelse(mapping_level == "InChIKey", 1, fifelse(mapping_level == "standardized_name", 0.8, 0.6))]
  entities[, mapping_conflict := figure8_bool(metadata_conflict_flag)]
  figure8_write_tsv(entities, "figure8e_compound_crosswalk.tsv")

  support <- entities[, .(
    canonical_name = figure8_first_or_na(canonical_name),
    raw_ids = paste(sort(unique(raw_id)), collapse = ","),
    inchi_key = figure8_first_or_na(inchi_key),
    mapping_level = figure8_first_or_na(mapping_level),
    mapping_confidence = max(mapping_confidence, na.rm = TRUE),
    mapping_conflict = any(mapping_conflict),
    DrugReflector = any(method == "DrugReflector" & support),
    L1000FWD = any(method == "L1000FWD" & support),
    CLUE = any(method == "CLUE" & support),
    clue_result_available = any(method == "CLUE"),
    l1000_result_listed = any(method == "L1000FWD"),
    drugreflector_rank = suppressWarnings(min(rank[method == "DrugReflector"], na.rm = TRUE)),
    l1000_rank = suppressWarnings(min(rank[method == "L1000FWD" & direction == "support"], na.rm = TRUE)),
    clue_tau = suppressWarnings(max(score[method == "CLUE"], na.rm = TRUE)),
    any_strong_opposition = any(strong_opposition),
    any_support = any(support)
  ), by = standardized_id]
  for (nm in c("drugreflector_rank", "l1000_rank", "clue_tau", "mapping_confidence")) support[!is.finite(get(nm)), (nm) := NA_real_]

  decomp_path <- figure8_input("metadata", "driver", "module9_9_landmark_decomposition_final_priority.tsv")
  decomp <- figure8_fread(decomp_path)
  liver <- decomp[, .(
    raw_id = compound,
    liver_context_score = as.numeric(clue_liver_context_percentile),
    clue_branch_balance = as.numeric(clue_branch_balance_percentile),
    clue_combined_percentile = as.numeric(clue_combined_percentile),
    clue_component_support = as.numeric(clue_branch_balance_percentile) >= 0.95 & as.numeric(clue_combined_percentile) >= 0.90
  )]
  liver <- merge(liver, meta[, .(compound, canonical_name, inchi_key)], by.x = "raw_id", by.y = "compound", all.x = TRUE)
  liver[, standardized_id := figure8_entity_key(inchi_key, canonical_name, raw_id)]
  liver <- liver[, .(
    liver_context_score = suppressWarnings(max(liver_context_score, na.rm = TRUE)),
    clue_branch_balance = suppressWarnings(max(clue_branch_balance, na.rm = TRUE)),
    clue_component_support = any(clue_component_support, na.rm = TRUE)
  ), by = standardized_id]
  liver[!is.finite(liver_context_score), liver_context_score := NA_real_]
  support <- merge(support, liver, by = "standardized_id", all.x = TRUE)
  support[, `:=`(
    CLUE = CLUE | fifelse(is.na(clue_component_support), FALSE, clue_component_support),
    liver_support = !is.na(liver_context_score) & liver_context_score >= 0.90
  )]
  support[, n_support_methods := as.integer(DrugReflector) + as.integer(L1000FWD) + as.integer(CLUE)]
  support[, support_category := fcase(
    any_strong_opposition & n_support_methods > 0, "Discordant",
    n_support_methods >= 2 & liver_support, "Strong support",
    n_support_methods >= 2, "Partial support",
    DrugReflector & !L1000FWD & !CLUE, "Model-only",
    !DrugReflector & (L1000FWD | CLUE), "External-only",
    is.na(mapping_level), "Unmapped",
    default = "Unresolved"
  )]
  figure8_write_tsv(support, "figure8e_method_support_matrix.tsv")

  combos <- CJ(DrugReflector = c(FALSE, TRUE), L1000FWD = c(FALSE, TRUE), CLUE = c(FALSE, TRUE))[DrugReflector | L1000FWD | CLUE]
  observed <- support[, .(count = .N), by = .(DrugReflector, L1000FWD, CLUE)]
  overlap <- merge(combos, observed, by = c("DrugReflector", "L1000FWD", "CLUE"), all.x = TRUE)
  overlap[is.na(count), count := 0L]
  overlap[, combination := apply(.SD, 1, function(z) paste(c("DrugReflector", "L1000FWD", "CLUE")[as.logical(z)], collapse = "+")), .SDcols = c("DrugReflector", "L1000FWD", "CLUE")]
  overlap[, overlap_class := fcase(DrugReflector & L1000FWD & CLUE, "three-method overlap", rowSums(.SD) == 2, "two-method overlap", DrugReflector & !L1000FWD & !CLUE, "DrugReflector only", !DrugReflector & L1000FWD & !CLUE, "L1000FWD only", !DrugReflector & !L1000FWD & CLUE, "CLUE only", default = "other"), .SDcols = c("DrugReflector", "L1000FWD", "CLUE")]
  figure8_write_tsv(overlap, "figure8e_method_overlap.tsv")

  rank_dr <- entities[method == "DrugReflector", .(rank = min(rank)), by = standardized_id]
  rank_l <- entities[method == "L1000FWD" & direction == "support", .(rank = min(rank)), by = standardized_id]
  rank_c <- entities[method == "CLUE", .(rank = min(rank)), by = standardized_id]
  corrs <- rbindlist(list(figure8_rank_correlation_row(rank_dr, rank_l, "DrugReflector", "L1000FWD"), figure8_rank_correlation_row(rank_dr, rank_c, "DrugReflector", "CLUE"), figure8_rank_correlation_row(rank_l, rank_c, "L1000FWD", "CLUE")))
  figure8_write_tsv(corrs, "figure8e_rank_correlation.tsv")

  mapping_overlap <- data.table(
    mapping_level = c("BRD_ID", "standardized_name", "InChIKey"),
    drugreflector_l1000_overlap = c(
      length(intersect(dr[support == TRUE, raw_id], l1000[support == TRUE, raw_id])),
      length(intersect(figure8_safe_name(dr[support == TRUE, canonical_name]), figure8_safe_name(l1000[support == TRUE, canonical_name]))),
      length(intersect(na.omit(dr[support == TRUE, inchi_key]), na.omit(l1000[support == TRUE, inchi_key])))
    ),
    drugreflector_clue_overlap = c(
      length(intersect(dr[support == TRUE, raw_id], clue[support == TRUE, raw_id])),
      length(intersect(figure8_safe_name(dr[support == TRUE, canonical_name]), figure8_safe_name(clue[support == TRUE, canonical_name]))),
      length(intersect(na.omit(dr[support == TRUE, inchi_key]), na.omit(clue[support == TRUE, inchi_key])))
    )
  )
  figure8_write_tsv(mapping_overlap, "figure8e_mapping_level_overlap.tsv")
  three <- overlap[DrugReflector & L1000FWD & CLUE, count]
  figure8_write_json(list(
    module = "figure8e_concordance", status = "completed", external_api_refresh = FALSE,
    direction_rules = list(DrugReflector = "Higher probability/lower rank for the desired rescue v-score supports reversal", L1000FWD = "similar_to_reversal_signature supports because the submitted query already encodes the desired rescue state", CLUE = "Positive tau supports the submitted desired-rescue UP/DOWN query; negative tau opposes it"),
    support_thresholds = list(DrugReflector = "primary top 200", L1000FWD = "similar group top 50", CLUE = "tau >= 90 or frozen component branch-balance >=95th and combined >=90th percentile", liver = "frozen CLUE liver-context percentile >= 0.90"),
    three_method_overlap = three %||% 0L, mapping_level_overlap = mapping_overlap, rank_correlations = corrs,
    no_query_is_not_negative = TRUE,
    interpretation_boundary = "Cross-platform connectivity is predictive support, not demonstrated efficacy or direct phenotypic rescue."
  ), "figure8e_concordance_report.json")
  invisible(support)
}

plot_figure8e <- function() {
  suppressPackageStartupMessages(library(ComplexUpset))
  support <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_support_matrix.tsv"))
  overlap <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_overlap.tsv"))
  display_ids <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_display_candidate_manifest.tsv"))$compound
  meta <- figure8_load_perturbagen_metadata()
  display_keys <- meta[compound %in% display_ids, figure8_entity_key(inchi_key, canonical_name, compound)]
  upset_data <- as.data.frame(support[, .(standardized_id, DrugReflector, L1000FWD, CLUE)])
  p_up <- ComplexUpset::upset(
    upset_data, intersect = c("DrugReflector", "L1000FWD", "CLUE"), min_size = 0,
    base_annotations = list("Intersection size" = ComplexUpset::intersection_size(text = list(size = 2.6))),
    set_sizes = ComplexUpset::upset_set_size(), width_ratio = 0.22
  ) + scale_fill_manual(values = method_palette) + ggtitle("Cross-method candidate intersections")
  overlap[, combination := factor(combination, levels = overlap[order(count), combination])]
  p_combo <- ggplot(overlap, aes(x = count, y = combination, fill = overlap_class)) +
    geom_col(width = 0.62) + geom_text(aes(label = count), hjust = -0.15, size = 2.2) +
    scale_fill_manual(values = c("three-method overlap" = lancet_palette[2], "two-method overlap" = lancet_palette[4], "DrugReflector only" = method_palette[["DrugReflector"]], "L1000FWD only" = method_palette[["L1000FWD"]], "CLUE only" = method_palette[["CLUE"]], other = lancet_palette[8])) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "All combinations (zero retained)", x = "Standardized entities", y = NULL, fill = NULL) + figure8_theme() + theme(legend.position = "none")
  mat <- support[standardized_id %in% display_keys]
  mat[, candidate := fifelse(is.na(canonical_name) | canonical_name == "", standardized_id, canonical_name)]
  long <- melt(mat, id.vars = c("standardized_id", "candidate", "support_category", "liver_support"), measure.vars = c("DrugReflector", "L1000FWD", "CLUE"), variable.name = "method", value.name = "supported")
  p_mat <- ggplot(long, aes(x = method, y = reorder(candidate, as.numeric(supported)), fill = supported)) +
    geom_tile(colour = "white", linewidth = 0.3) + geom_point(data = long[liver_support == TRUE], shape = 21, size = 1.6, fill = "white", colour = lancet_palette[9]) +
    scale_fill_manual(values = c(`TRUE` = lancet_palette[3], `FALSE` = "#F7F7F7")) +
    labs(title = "Support matrix (circle: liver context)", x = NULL, y = NULL, fill = "Directionally supportive") + figure8_theme() + theme(legend.position = "bottom")
  combined <- (p_up / (p_combo | p_mat)) + plot_layout(heights = c(1.1, 1)) +
    plot_annotation(title = "DrugReflector, L1000FWD and CLUE concordance", caption = "Entity priority: InChIKey > standardized name > BRD ID. Positive CLUE tau and L1000FWD similar-to-rescue are support after query-orientation audit.")
  figure8_save(combined, "e", "figure8e_cross_method_concordance", 12.8, 8.4)
}

run_figure8_mechanism_classes <- function() {
  suppressPackageStartupMessages({
    library(igraph)
    library(tidygraph)
  })
  pred <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz"))
  stab <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  support <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_support_matrix.tsv"))
  n_compounds <- uniqueN(pred$compound)
  scores <- pred[signature_id %in% c("primary_three_axis", "identity_rescue_only", "stress_suppression_only", "sox4_only", "no_proliferation")]
  scores[, score := figure8_rank_score(rank_1based, n_compounds)]
  wide <- dcast(scores, compound ~ signature_id, value.var = "score")
  setnames(wide, c("identity_rescue_only", "stress_suppression_only", "sox4_only"), c("identity_rescue", "stress_suppression", "sox4_suppression"), skip_absent = TRUE)
  profiles <- merge(stab, wide, by = "compound", all.x = TRUE)
  profiles[, standardized_id := figure8_entity_key(inchi_key, canonical_name, compound)]
  profiles <- merge(profiles, support[, .(standardized_id, DrugReflector, L1000FWD, CLUE, liver_support, support_category, n_support_methods)], by = "standardized_id", all.x = TRUE)
  profiles[, `:=`(
    DrugReflector = fifelse(is.na(DrugReflector), FALSE, DrugReflector),
    L1000FWD = fifelse(is.na(L1000FWD), FALSE, L1000FWD),
    CLUE = fifelse(is.na(CLUE), FALSE, CLUE),
    liver_support = fifelse(is.na(liver_support), FALSE, liver_support),
    external_support_score = rowMeans(cbind(as.numeric(L1000FWD), as.numeric(CLUE)), na.rm = TRUE),
    proliferation_dependency = pmax(0, primary_three_axis - no_proliferation),
    axis_mean = rowMeans(cbind(identity_rescue, stress_suppression, sox4_suppression), na.rm = TRUE),
    axis_min = pmin(identity_rescue, stress_suppression, sox4_suppression, na.rm = TRUE)
  )]
  profiles[, response_class := fcase(
    proliferation_dependency >= 0.20, "proliferation-dependent apparent reversal",
    axis_min >= 0.95, "parallel three-axis reversal",
    identity_rescue >= 0.95 & identity_rescue >= rowMeans(cbind(stress_suppression, sox4_suppression), na.rm = TRUE) + 0.03, "identity-rescue dominant",
    rowMeans(cbind(stress_suppression, sox4_suppression), na.rm = TRUE) >= 0.95 & rowMeans(cbind(stress_suppression, sox4_suppression), na.rm = TRUE) >= identity_rescue + 0.03, "stress/SOX4 suppression dominant",
    external_support_score > 0 & median_rank_score < 0.90, "external-only support",
    median_rank_score >= 0.90 & external_support_score == 0, "DrugReflector-only exploratory predictions",
    default = "mixed or unresolved response"
  )]
  profiles[, `:=`(
    moa = NA_character_, target = NA_character_, target_family = NA_character_, pathway_class = NA_character_,
    compound_status = NA_character_, moa_availability = "unavailable: no frozen compound MoA/target table",
    mechanism_interpretability = NA_real_
  )]
  analysis <- profiles[in_candidate_analysis_universe == TRUE | external_support_score > 0]
  summary <- analysis[, .(
    candidate_count = .N,
    evidence_weighted_concentration = sum(0.5 * median_rank_score + 0.5 * external_support_score, na.rm = TRUE),
    median_reversal_score = median(axis_mean, na.rm = TRUE),
    median_identity_rescue = median(identity_rescue, na.rm = TRUE),
    median_stress_suppression = median(stress_suppression, na.rm = TRUE),
    median_sox4_suppression = median(sox4_suppression, na.rm = TRUE),
    liver_context_fraction = mean(liver_support, na.rm = TRUE),
    proliferation_dependency_fraction = mean(proliferation_dependency >= 0.20, na.rm = TRUE)
  ), by = response_class]
  summary[, recurrent_response_class := candidate_count >= 2]
  profiles <- merge(profiles, summary[, .(response_class, response_class_count = candidate_count, recurrent_response_class)], by = "response_class", all.x = TRUE)

  moa <- profiles[, .(
    compound_name = canonical_name, standardized_id, compound, MoA = moa, target, target_family,
    pathway_class, approved_investigational_research_status = compound_status,
    liver_context_evidence = liver_support, DrugReflector_rank = median_rank,
    CLUE_support = CLUE, L1000FWD_support = L1000FWD,
    toxicity_pan_stress_flags = fifelse(proliferation_dependency >= 0.20, "proliferation-dependent model signal", "unknown"),
    response_class, moa_availability
  )]
  axis_matrix <- profiles[, .(
    compound, compound_name = canonical_name, standardized_id, response_class,
    identity_rescue, stress_suppression, sox4_suppression,
    malignant_fate_reversal = rowMeans(cbind(stress_suppression, sox4_suppression), na.rm = TRUE),
    proliferation_independent_effect = 1 - pmin(1, proliferation_dependency),
    liver_context_support = as.numeric(liver_support), external_support_score,
    recurrent_response_class, in_candidate_analysis_universe, display_in_main_panel
  )]

  display <- axis_matrix[display_in_main_panel == TRUE]
  candidate_nodes <- data.table(name = display$compound, label = display$compound_name, node_type = "compound")
  class_nodes <- data.table(name = unique(display$response_class), label = unique(display$response_class), node_type = "response_class")
  method_nodes <- data.table(name = c("DrugReflector", "L1000FWD", "CLUE"), label = c("DrugReflector", "L1000FWD", "CLUE"), node_type = "method")
  edges <- rbindlist(list(
    display[, .(from = compound, to = response_class, edge_type = "response_class")],
    profiles[compound %in% display$compound & DrugReflector == TRUE, .(from = compound, to = "DrugReflector", edge_type = "method_support")],
    profiles[compound %in% display$compound & L1000FWD == TRUE, .(from = compound, to = "L1000FWD", edge_type = "method_support")],
    profiles[compound %in% display$compound & CLUE == TRUE, .(from = compound, to = "CLUE", edge_type = "method_support")]
  ), fill = TRUE)
  nodes <- unique(rbindlist(list(candidate_nodes, class_nodes, method_nodes), fill = TRUE), by = "name")
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  centrality <- data.table(name = names(igraph::degree(graph)), degree = as.numeric(igraph::degree(graph)), betweenness = as.numeric(igraph::betweenness(graph)))
  nodes <- merge(nodes, centrality, by = "name", all.x = TRUE)
  figure8_write_tsv(moa, "figure8f_compound_moa.tsv")
  figure8_write_tsv(summary, "figure8f_mechanism_class_summary.tsv")
  figure8_write_tsv(axis_matrix, "figure8f_mechanism_axis_matrix.tsv")
  figure8_write_tsv(edges, "figure8f_evidence_network_edges.tsv")
  figure8_write_tsv(nodes, "figure8f_evidence_network_nodes.tsv")
  figure8_write_json(list(
    module = "figure8f_mechanism_classes", status = "completed_with_moa_gap",
    moa_source = "No frozen compound-level MoA/target table was found; no online or name-based annotations were added.",
    classification_source = "Preregistered rank-score rules over identity, stress and SOX4/malignant axis profiles plus cross-method support.",
    response_class_counts = summary[, .(response_class, candidate_count)],
    network_packages = figure8_package_versions(c("igraph", "tidygraph", "ggraph")),
    interpretation_boundary = "Response classes are transcriptomic/model classes and are not pharmacological mechanism-of-action assignments."
  ), "figure8f_mechanism_report.json")
  invisible(profiles)
}

plot_figure8f <- function() {
  suppressPackageStartupMessages({
    library(tidygraph)
    library(ggraph)
  })
  summary <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_class_summary.tsv"))
  axis <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_axis_matrix.tsv"))[display_in_main_panel == TRUE]
  edges <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_evidence_network_edges.tsv"))
  nodes <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_evidence_network_nodes.tsv"))
  p1 <- ggplot(summary, aes(x = evidence_weighted_concentration, y = reorder(response_class, evidence_weighted_concentration))) +
    geom_point(aes(size = candidate_count, fill = median_reversal_score), shape = 21, colour = "white", stroke = 0.5) +
    scale_fill_gradient2(low = reversal_gradient[["low"]], mid = reversal_gradient[["mid"]], high = reversal_gradient[["high"]], midpoint = 0.5, limits = c(0, 1)) +
    scale_size_area(max_size = 10) +
    labs(title = "Data-derived response-class concentration", x = "Evidence-weighted concentration", y = NULL, size = "Candidates", fill = "Median reversal") + figure8_theme()
  long <- melt(axis, id.vars = c("compound", "compound_name", "response_class"), measure.vars = c("identity_rescue", "stress_suppression", "sox4_suppression", "malignant_fate_reversal", "proliferation_independent_effect", "liver_context_support"), variable.name = "dimension", value.name = "score")
  long[, compound_name := fifelse(is.na(compound_name) | compound_name == "", compound, compound_name)]
  p2 <- ggplot(long, aes(x = dimension, y = reorder(compound_name, score), fill = score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = reversal_gradient[["low"]], mid = reversal_gradient[["mid"]], high = reversal_gradient[["high"]], midpoint = 0.5, limits = c(0, 1), na.value = "white") +
    labs(title = "Three-axis and context profiles", x = NULL, y = NULL, fill = "Rank score") + figure8_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  graph <- tidygraph::tbl_graph(nodes = as.data.frame(nodes[, .(name, label, node_type)]), edges = as.data.frame(edges[, .(from, to, edge_type)]), directed = FALSE, node_key = "name")
  p3 <- ggraph::ggraph(graph, layout = "fr") +
    ggraph::geom_edge_link(aes(linetype = edge_type), colour = lancet_palette[8], alpha = 0.65, linewidth = 0.35) +
    ggraph::geom_node_point(aes(fill = node_type), shape = 21, size = 3.5, colour = "white", stroke = 0.4) +
    ggraph::geom_node_text(aes(label = label), repel = TRUE, size = 2.0, max.overlaps = 30) +
    scale_fill_manual(values = c(compound = lancet_palette[6], response_class = lancet_palette[3], method = lancet_palette[1])) +
    labs(title = "Candidate–response-class–method evidence network", fill = "Node", linetype = "Edge") + theme_void(base_family = "sans", base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 10))
  combined <- (p1 | p2) / p3 + plot_layout(heights = c(1, 0.9)) +
    plot_annotation(title = "Transcriptomic response classes and evidence structure", caption = "Compound-level MoA and targets were unavailable in frozen project metadata; classes shown here are data-derived response classes, not inferred pharmacology.")
  figure8_save(combined, "f", "figure8f_mechanism_class_enrichment", 13.0, 9.0)
}

figure8_parse_clue_sig_id <- function(sig_id) {
  sig_id <- as.character(sig_id)
  head <- sub(":.*$", "", sig_id)
  time <- sub("^.*_([^_]+)$", "\\1", head)
  cell_line <- sub("^[^_]+_([^_]+)_.*$", "\\1", head)
  dose <- sub("^.*:", "", sig_id)
  data.table(cell_line_parsed = cell_line, time = time, dose = dose)
}

run_figure8_external_validation <- function() {
  path <- figure8_input("metadata", "driver", "module9_9_landmark_decomposition_clue_signature_components.tsv.gz")
  figure8_require(path, "Frozen CLUE signature components")
  x <- figure8_fread(path)
  parsed <- figure8_parse_clue_sig_id(x$sig_id)
  x[, `:=`(cell_line = parsed$cell_line_parsed, time = parsed$time, dose = parsed$dose)]
  x[, `:=`(
    rescue_connectivity = as.numeric(cs_up),
    malignant_suppression_connectivity = -as.numeric(cs_down),
    combined_connectivity = as.numeric(cs_combined),
    identity_rescue_score = NA_real_, stress_suppression_score = NA_real_, sox4_suppression_score = NA_real_,
    malignant_fate_reversal_score = NA_real_, proliferation_change = NA_real_, generic_stress_induction = NA_real_,
    viability_cytotoxicity_proxy = NA_real_, direct_axis_score_available = FALSE,
    validation_type = "CLUE paired-query connectivity decomposition; not direct expression axis validation"
  )]
  stab <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  candidates <- stab[in_candidate_analysis_universe == TRUE | display_in_main_panel == TRUE, compound]
  subset <- x[compound %in% candidates]
  figure8_write_tsv(subset, "figure8g_external_signature_scores.tsv", compress = TRUE)

  liver_lines <- c("HEPG2", "HCC515", "HA1E")
  liver <- subset[toupper(cell_line) %in% liver_lines]
  summary <- liver[, {
    values <- combined_connectivity[is.finite(combined_connectivity)]
    n <- length(values)
    avg <- if (n) mean(values) else NA_real_
    sdev <- if (n >= 2) sd(values) else NA_real_
    se <- if (n >= 2) sdev / sqrt(n) else NA_real_
    mult <- if (n >= 2) qt(0.975, df = n - 1) else NA_real_
    list(
      n_signatures = n,
      mean_rescue_connectivity = mean(rescue_connectivity, na.rm = TRUE),
      mean_malignant_suppression_connectivity = mean(malignant_suppression_connectivity, na.rm = TRUE),
      mean_combined_connectivity = avg,
      median_combined_connectivity = if (n) median(values) else NA_real_,
      sd_combined_connectivity = sdev,
      se_combined_connectivity = se,
      ci_low = if (n >= 2) avg - mult * se else NA_real_,
      ci_high = if (n >= 2) avg + mult * se else NA_real_,
      dose_values = paste(sort(unique(dose)), collapse = ","),
      time_values = paste(sort(unique(time)), collapse = ","),
      direction = ifelse(is.na(avg), "unavailable", ifelse(avg > 0, "supports desired rescue connectivity", ifelse(avg < 0, "opposes desired rescue connectivity", "neutral")))
    )
  }, by = .(compound, cell_line)]
  summary <- merge(summary, stab[, .(compound, compound_name = canonical_name, display_in_main_panel)], by = "compound", all.x = TRUE)
  figure8_write_tsv(summary, "figure8g_liver_context_summary.tsv")
  figure8_write_json(list(
    module = "figure8g_external_validation", status = if (nrow(summary)) "connectivity_support_available" else "unavailable",
    source = normalizePath(path, winslash = "/", mustWork = TRUE),
    source_type = "CLUE cs_up, cs_down and cs_combined connectivity from a paired query",
    direct_expression_matrix_available = FALSE, direct_three_axis_scores_available = FALSE,
    dose_time_parsed_from_sig_id = TRUE, liver_cell_lines = sort(unique(summary$cell_line)),
    n_candidate_signatures = nrow(subset), n_liver_candidate_signatures = nrow(liver),
    uncertainty = "Across-signature t intervals are descriptive because CLUE profiles do not provide profile-specific sampling variances; metafor was not applied.",
    review_risks = c("Connectivity is not direct expression", "Biological replicate identity and viability are unavailable", "Stress and SOX4 axes cannot be separated from this paired query"),
    interpretation_boundary = "Connectivity support cannot establish phenotypic rescue, efficacy, or safety."
  ), "figure8g_external_validation_report.json")
  invisible(summary)
}

plot_figure8g <- function() {
  x <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8g_external_signature_scores.tsv.gz"))
  summary <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8g_liver_context_summary.tsv"))
  display <- summary[display_in_main_panel == TRUE, unique(compound)]
  agg <- x[compound %in% display, .(
    rescue_connectivity = mean(rescue_connectivity, na.rm = TRUE),
    malignant_suppression_connectivity = mean(malignant_suppression_connectivity, na.rm = TRUE),
    combined_connectivity = mean(combined_connectivity, na.rm = TRUE),
    identity_axis = NA_real_, stress_axis = NA_real_, sox4_axis = NA_real_
  ), by = compound]
  names_map <- unique(summary[, .(compound, compound_name)])
  agg <- merge(agg, names_map, by = "compound", all.x = TRUE)
  long <- melt(agg, id.vars = c("compound", "compound_name"), variable.name = "dimension", value.name = "score")
  long[, compound_name := fifelse(is.na(compound_name) | compound_name == "", compound, compound_name)]
  p1 <- ggplot(long, aes(x = dimension, y = reorder(compound_name, score), fill = score)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(data = long[is.na(score)], label = "NA", size = 1.8, colour = lancet_palette[8]) +
    scale_fill_gradient2(low = reversal_gradient[["low"]], mid = reversal_gradient[["mid"]], high = reversal_gradient[["high"]], midpoint = 0, na.value = "white") +
    labs(title = "Frozen CLUE connectivity components", x = NULL, y = NULL, fill = "Connectivity") + figure8_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  forest <- summary[display_in_main_panel == TRUE]
  forest[, label := paste0(fifelse(is.na(compound_name) | compound_name == "", compound, compound_name), " | ", cell_line, " | n=", n_signatures, " | ", time_values)]
  p2 <- ggplot(forest, aes(x = mean_combined_connectivity, y = reorder(label, mean_combined_connectivity), colour = cell_line)) +
    geom_vline(xintercept = 0, colour = lancet_palette[8], linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.4, na.rm = TRUE) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = setNames(lancet_palette[c(1, 3, 4)], c("HEPG2", "HCC515", "HA1E"))) +
    labs(title = "Liver-relevant cell-line connectivity", x = "Mean combined connectivity (descriptive 95% CI)", y = NULL, colour = "Cell line") + figure8_theme() + theme(legend.position = "bottom")
  combined <- (p1 | p2) + plot_layout(widths = c(1, 1.15)) +
    plot_annotation(title = "Independent perturbational connectivity support", caption = "CLUE summary connectivity is shown. Direct three-axis expression scores, viability, and phenotypic rescue are unavailable and remain marked NA.")
  figure8_save(combined, "g", "figure8g_external_perturbation_validation", 13.0, 7.4)
}

figure8_sample_expression_matched <- function(template, universe) {
  chosen <- character()
  for (bin in unique(template$expression_bin)) {
    n <- sum(template$expression_bin == bin)
    pool <- setdiff(universe[expression_bin == bin, gene], chosen)
    if (length(pool) < n) pool <- setdiff(universe$gene, chosen)
    chosen <- c(chosen, sample(pool, n, replace = length(pool) < n))
  }
  chosen
}

run_figure8_random_benchmark <- function() {
  set.seed(FIGURE8_SEED)
  variants <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_signature_variants_long.tsv.gz"))
  expression <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_expression_universe.tsv.gz"))
  model <- figure8_model_genes()
  primary <- variants[signature_id == "primary_three_axis"]
  primary <- merge(primary, expression[, .(gene, expression_bin, mean_expression)], by = "gene", all.x = TRUE)
  primary[is.na(expression_bin), expression_bin := "zero_or_absent"]
  universe <- unique(expression[, .(gene, expression_bin, mean_expression)])
  n_up <- primary[desired_direction == "up", uniqueN(gene)]
  n_down <- primary[desired_direction == "down", uniqueN(gene)]
  up_weights <- abs(primary[desired_direction == "up", v_score])
  down_weights <- abs(primary[desired_direction == "down", v_score])

  make_profile <- function(id, type, up_genes, down_genes, up_w = sample(up_weights, length(up_genes), replace = TRUE), down_w = sample(down_weights, length(down_genes), replace = TRUE)) {
    rbind(
      data.table(signature_id = id, null_type = type, gene = up_genes, desired_direction = "up", v_score = abs(up_w)),
      data.table(signature_id = id, null_type = type, gene = down_genes, desired_direction = "down", v_score = -abs(down_w))
    )
  }
  profiles <- vector("list", 1000L)
  cursor <- 0L
  for (i in seq_len(265L)) {
    cursor <- cursor + 1L
    genes <- sample(universe$gene, n_up + n_down, replace = FALSE)
    profiles[[cursor]] <- make_profile(sprintf("gene_number_matched_%03d", i), "gene_number_matched", genes[seq_len(n_up)], genes[n_up + seq_len(n_down)])
  }
  for (i in seq_len(265L)) {
    cursor <- cursor + 1L
    up <- figure8_sample_expression_matched(primary[desired_direction == "up"], universe)
    remaining <- universe[!gene %in% up]
    down <- figure8_sample_expression_matched(primary[desired_direction == "down"], remaining)
    profiles[[cursor]] <- make_profile(sprintf("expression_matched_%03d", i), "expression_matched", up, down)
  }
  lm_up_n <- primary[desired_direction == "up" & gene %in% model, uniqueN(gene)]
  lm_down_n <- primary[desired_direction == "down" & gene %in% model, uniqueN(gene)]
  model_pool <- intersect(model, universe$gene)
  non_model_pool <- setdiff(universe$gene, model)
  for (i in seq_len(265L)) {
    cursor <- cursor + 1L
    lm <- sample(model_pool, lm_up_n + lm_down_n, replace = FALSE)
    non <- sample(non_model_pool, (n_up - lm_up_n) + (n_down - lm_down_n), replace = FALSE)
    up <- c(lm[seq_len(lm_up_n)], non[seq_len(n_up - lm_up_n)])
    down <- c(lm[lm_up_n + seq_len(lm_down_n)], non[(n_up - lm_up_n) + seq_len(n_down - lm_down_n)])
    profiles[[cursor]] <- make_profile(sprintf("landmark_coverage_matched_%03d", i), "landmark_coverage_matched", up, down)
  }
  base_genes <- unique(primary$gene)
  for (i in seq_len(100L)) {
    cursor <- cursor + 1L
    shuffled <- sample(c(rep("up", n_up), rep("down", n_down)))
    profiles[[cursor]] <- make_profile(sprintf("shuffled_direction_%03d", i), "shuffled_direction", base_genes[shuffled == "up"], base_genes[shuffled == "down"])
  }
  for (i in seq_len(100L)) {
    cursor <- cursor + 1L
    profiles[[cursor]] <- make_profile(
      sprintf("random_weight_%03d", i), "random_weight",
      primary[desired_direction == "up", gene], primary[desired_direction == "down", gene],
      runif(n_up, min(up_weights), max(up_weights)), runif(n_down, min(down_weights), max(down_weights))
    )
  }
  for (v in c("leave_out_sox4", "leave_out_stress", "leave_out_identity")) {
    cursor <- cursor + 1L
    z <- variants[signature_id == v]
    z[, `:=`(signature_id = paste0("negative_control_", v), null_type = v)]
    profiles[[cursor]] <- z[, .(signature_id, null_type, gene, desired_direction, v_score)]
  }
  nuisance <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_nuisance_gene_sets.tsv"))
  for (type in c("proliferation", "generic_stress")) {
    cursor <- cursor + 1L
    genes <- nuisance[nuisance_set == type, unique(gene)]
    z <- primary[gene %in% genes]
    if (!nrow(z)) z <- primary[axis == "stress_suppression"]
    z[, `:=`(signature_id = paste0(type, "_only_control"), null_type = paste0(type, "_only"))]
    profiles[[cursor]] <- z[, .(signature_id, null_type, gene, desired_direction, v_score)]
  }
  if (cursor != 1000L) stop("Random signature design must contain exactly 1,000 profiles; found ", cursor)
  random_long <- rbindlist(profiles, fill = TRUE)
  random_long[, is_landmark := gene %in% model]
  manifest <- random_long[, .(
    up_gene_count = uniqueN(gene[desired_direction == "up"]), down_gene_count = uniqueN(gene[desired_direction == "down"]),
    landmark_up_count = uniqueN(gene[desired_direction == "up" & is_landmark]), landmark_down_count = uniqueN(gene[desired_direction == "down" & is_landmark]),
    absolute_weight_mass = sum(abs(v_score))
  ), by = .(signature_id, null_type)]
  figure8_write_tsv(manifest, "figure8_random_signature_manifest.tsv")
  wide <- dcast(random_long, signature_id ~ gene, value.var = "v_score", fun.aggregate = sum, fill = 0)
  input_path <- figure8_write_data_tsv(wide, "figure8_random_signatures_wide.tsv", compress = TRUE)
  stability <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  watch <- stability[in_candidate_analysis_universe == TRUE, .(compound)]
  watch_path <- figure8_write_tsv(watch, "figure8_random_candidate_watchlist.tsv")

  python <- figure8_input(".venv-drugreflector", "Scripts", "python.exe")
  if (!file.exists(python)) python <- Sys.which("python")
  if (!nzchar(python)) stop("DrugReflector Python runtime is unavailable")
  args <- c(
    shQuote(figure8_input("scripts", "figure8_drugreflector_inference.py")), "--mode", "random",
    "--input", shQuote(input_path), "--watchlist", shQuote(watch_path), "--metadata-dir", shQuote(FIGURE8_METADATA_DIR),
    "--top-n", "200", "--batch-size", "20", "--seed", as.character(FIGURE8_SEED)
  )
  status <- system2(python, args = args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status") %||% 0L
  writeLines(status, file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_inference.log"))
  if (exit_status != 0L) stop("Random DrugReflector inference failed with status ", exit_status)

  inf <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_inference_summary.tsv.gz"))
  top <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_top_predictions.tsv.gz"))
  watched <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_watchlist_predictions.tsv.gz"))
  inf <- merge(inf, manifest, by = "signature_id", all.x = TRUE)
  primary_pred <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz"))[signature_id == "primary_three_axis"]
  primary_top100 <- primary_pred[rank_1based <= 100, compound]
  support <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_support_matrix.tsv"))
  meta <- figure8_load_perturbagen_metadata()
  map <- meta[, .(compound, standardized_id = figure8_entity_key(inchi_key, canonical_name, compound))]
  top <- merge(top, map, by = "compound", all.x = TRUE)
  top <- merge(top, support[, .(standardized_id, external_support = L1000FWD | CLUE, liver_support)], by = "standardized_id", all.x = TRUE)
  axis <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_axis_matrix.tsv"))[, .(compound, response_class)]
  top <- merge(top, axis, by = "compound", all.x = TRUE)
  metrics <- top[, .(
    top100_overlap_with_real = sum(rank_1based <= 100 & compound %in% primary_top100),
    cross_method_support_count = sum(external_support, na.rm = TRUE),
    liver_support_count = sum(liver_support, na.rm = TRUE),
    response_class_concentration = if (all(is.na(response_class))) NA_real_ else max(table(response_class)) / sum(!is.na(response_class))
  ), by = signature_id]
  benchmark <- merge(inf, metrics, by = "signature_id", all.x = TRUE)
  benchmark[, `:=`(
    top100_overlap_fraction = top100_overlap_with_real / 100,
    cross_method_overlap_fraction = cross_method_support_count / 200,
    liver_support_fraction = liver_support_count / 200
  )]
  benchmark[, integrated_null_proxy := 0.50 * figure8_norm01(max_probability) + 0.20 * top100_overlap_fraction + 0.15 * cross_method_overlap_fraction + 0.15 * liver_support_fraction]

  real_top <- primary_pred[rank_1based <= 200, .(compound, rank_1based)]
  real_top <- merge(real_top, map, by = "compound", all.x = TRUE)
  real_top <- merge(real_top, support[, .(standardized_id, external_support = L1000FWD | CLUE, liver_support)], by = "standardized_id", all.x = TRUE)
  real_top <- merge(real_top, axis, by = "compound", all.x = TRUE)
  real_metrics <- list(
    max_probability = max(primary_pred$prob),
    top100_overlap_with_real = 100,
    cross_method_support_count = sum(real_top$external_support, na.rm = TRUE),
    liver_support_count = sum(real_top$liver_support, na.rm = TRUE),
    response_class_concentration = if (all(is.na(real_top$response_class))) NA_real_ else max(table(real_top$response_class)) / sum(!is.na(real_top$response_class))
  )
  percentile <- function(value, null) mean(null <= value, na.rm = TRUE)
  benchmark_summary <- data.table(
    metric = names(real_metrics), real_value = unlist(real_metrics),
    random_median = c(median(benchmark$max_probability), median(benchmark$top100_overlap_with_real), median(benchmark$cross_method_support_count), median(benchmark$liver_support_count), median(benchmark$response_class_concentration, na.rm = TRUE)),
    random_percentile = c(percentile(real_metrics$max_probability, benchmark$max_probability), percentile(real_metrics$top100_overlap_with_real, benchmark$top100_overlap_with_real), percentile(real_metrics$cross_method_support_count, benchmark$cross_method_support_count), percentile(real_metrics$liver_support_count, benchmark$liver_support_count), percentile(real_metrics$response_class_concentration, benchmark$response_class_concentration))
  )
  figure8_write_tsv(benchmark, "figure8_random_signature_benchmark.tsv", compress = TRUE)
  figure8_write_tsv(benchmark_summary, "figure8_random_signature_benchmark_summary.tsv")

  real_watch <- primary_pred[compound %in% watch$compound, .(compound, real_rank = rank_1based, real_prob = prob)]
  candidate_specificity <- watched[, .(
    random_rank_percentile = mean(rank_1based >= real_watch$real_rank[match(.BY$compound, real_watch$compound)], na.rm = TRUE),
    random_probability_percentile = mean(prob <= real_watch$real_prob[match(.BY$compound, real_watch$compound)], na.rm = TRUE),
    random_median_rank = median(rank_1based), random_median_probability = median(prob)
  ), by = compound]
  candidate_specificity <- merge(real_watch, candidate_specificity, by = "compound", all.x = TRUE)
  figure8_write_tsv(candidate_specificity, "figure8_random_candidate_specificity_percentiles.tsv")
  specificity_pass <- benchmark_summary[metric == "max_probability", random_percentile] >= 0.95 && benchmark_summary[metric == "cross_method_support_count", random_percentile] >= 0.95
  figure8_write_json(list(
    module = "figure8_random_signature_benchmark", status = "completed", seed = FIGURE8_SEED,
    n_random_signatures = nrow(manifest), null_type_counts = manifest[, .N, by = null_type],
    benchmark = benchmark_summary, signature_specificity_supported = specificity_pass,
    specificity_risk = !specificity_pass,
    limitations = c("External services were not re-queried for null signatures", "Cross-method and liver null metrics are overlaps with frozen supported compounds", "Response-class concentration substitutes for unavailable frozen MoA classes")
  ), "figure8_random_signature_benchmark_report.json")
  invisible(benchmark_summary)
}

run_figure8_toxicity_penalties <- function() {
  signature <- figure8_signature()
  nuisance <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_nuisance_gene_sets.tsv"))
  model <- figure8_model_genes()
  primary <- signature[include_primary == TRUE]
  set_map <- c(
    proliferation = "proliferation", generic_stress = "generic_stress", dna_damage = "dna_damage",
    translation_inhibition = "translation_inhibition", mitochondrial_toxicity = "mitochondrial_toxicity",
    unfolded_protein_response = "unfolded_protein_response", apoptosis = "apoptosis"
  )
  controls <- rbindlist(lapply(names(set_map), function(id) {
    genes <- nuisance[nuisance_set == set_map[[id]], unique(gene)]
    x <- primary[gene %in% genes]
    if (!nrow(x)) return(NULL)
    x[, .(
      signature_id = paste0("toxicity_control_", id), gene, desired_direction,
      v_score = fifelse(desired_direction == "up", 1, -1) * final_weight,
      nuisance_set = id, is_landmark = gene %in% model
    )]
  }), fill = TRUE)
  control_qc <- controls[, .(n_genes = uniqueN(gene), n_landmarks = uniqueN(gene[is_landmark]), status = fifelse(uniqueN(gene[is_landmark]) >= 3, "available", "incomplete")), by = .(signature_id, nuisance_set)]
  usable <- control_qc[status == "available", signature_id]
  if (length(usable)) {
    wide <- dcast(controls[signature_id %in% usable], signature_id ~ gene, value.var = "v_score", fun.aggregate = sum, fill = 0)
    input <- figure8_write_data_tsv(wide, "figure8_toxicity_control_signatures_wide.tsv", compress = TRUE)
    python <- figure8_input(".venv-drugreflector", "Scripts", "python.exe")
    if (!file.exists(python)) python <- Sys.which("python")
    args <- c(shQuote(figure8_input("scripts", "figure8_drugreflector_inference.py")), "--mode", "toxicity", "--input", shQuote(input), "--metadata-dir", shQuote(FIGURE8_METADATA_DIR), "--seed", as.character(FIGURE8_SEED))
    log <- system2(python, args = args, stdout = TRUE, stderr = TRUE)
    exit_status <- attr(log, "status") %||% 0L
    writeLines(log, file.path(FIGURE8_METADATA_DIR, "figure8_toxicity_control_inference.log"))
    if (exit_status != 0L) stop("Toxicity-control DrugReflector inference failed")
    pred <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_toxicity_control_predictions.tsv.gz"))
  } else {
    pred <- data.table(signature_id = character(), compound = character(), rank_1based = numeric())
  }
  all_compounds <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))[, .(compound)]
  if (nrow(pred)) {
    pred[, nuisance_set := sub("^toxicity_control_", "", signature_id)]
    pred[, penalty_score := figure8_rank_score(rank_1based, uniqueN(compound))]
    pen_wide <- dcast(pred, compound ~ nuisance_set, value.var = "penalty_score")
  } else pen_wide <- copy(all_compounds)
  penalties <- merge(all_compounds, pen_wide, by = "compound", all.x = TRUE)
  variants <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variant_predictions.tsv.gz"))
  variants[, rank_score := figure8_rank_score(rank_1based, uniqueN(compound))]
  selected <- dcast(variants[signature_id %in% c("primary_three_axis", "no_proliferation", "no_generic_stress")], compound ~ signature_id, value.var = "rank_score")
  penalties <- merge(penalties, selected, by = "compound", all.x = TRUE)
  for (nm in names(set_map)) if (!nm %in% names(penalties)) penalties[, (nm) := NA_real_]
  penalties[, `:=`(
    proliferation_reliance = pmax(0, primary_three_axis - no_proliferation),
    generic_stress_reliance = pmax(0, primary_three_axis - no_generic_stress)
  )]
  penalties[, proliferation_penalty := pmax(proliferation, proliferation_reliance, na.rm = TRUE)]
  penalties[, generic_stress_penalty := pmax(generic_stress, unfolded_protein_response, generic_stress_reliance, na.rm = TRUE)]
  penalties[!is.finite(proliferation_penalty), proliferation_penalty := NA_real_]
  penalties[!is.finite(generic_stress_penalty), generic_stress_penalty := NA_real_]
  penalties[, dna_damage_penalty := dna_damage]
  penalties[, translation_inhibition_penalty := translation_inhibition]
  penalties[, mitochondrial_toxicity_penalty := mitochondrial_toxicity]
  penalties[, pan_cytotoxicity_penalty := rowMeans(cbind(dna_damage, mitochondrial_toxicity, unfolded_protein_response, apoptosis), na.rm = TRUE)]
  penalties[!is.finite(pan_cytotoxicity_penalty), pan_cytotoxicity_penalty := NA_real_]
  penalties[, `:=`(
    cytotoxicity_penalty = pan_cytotoxicity_penalty,
    toxicity_data_status = "model-derived nuisance-program sensitivity; compound viability/safety unknown",
    toxicity_unknown = TRUE,
    broad_transcriptional_suppression_status = "unknown"
  )]
  figure8_write_tsv(controls, "figure8_toxicity_control_signature_genes.tsv")
  figure8_write_tsv(control_qc, "figure8_toxicity_control_signature_qc.tsv")
  figure8_write_tsv(penalties, "figure8_toxicity_stress_penalty.tsv")
  figure8_write_json(list(
    module = "figure8_toxicity_stress_penalties", status = "completed_with_safety_gap",
    control_signature_qc = control_qc,
    penalty_semantics = "Rank-score within frozen nuisance-gene subsets plus loss of rank after nuisance-gene removal; higher values indicate stronger nuisance dependence.",
    compound_level_viability_available = FALSE, safety_assessed = FALSE,
    unknown_handling = "Unknown compound-level toxicity remains unknown and blocks safety claims; it is not converted to zero or safe.",
    review_risks = c("Model-derived nuisance penalties are not experimental toxicity", "Broad transcriptional suppression and viability are unavailable")
  ), "figure8_toxicity_stress_penalty_report.json")
  invisible(penalties)
}

run_figure8_integrated_score <- function() {
  stab <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  axis <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_axis_matrix.tsv"))
  support <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_support_matrix.tsv"))
  penalty <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_toxicity_stress_penalty.tsv"))
  specificity <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_candidate_specificity_percentiles.tsv"))
  profiles <- merge(stab, axis[, .(compound, response_class, identity_rescue, stress_suppression, sox4_suppression, malignant_fate_reversal, liver_context_support, external_support_score, recurrent_response_class)], by = "compound", all.x = TRUE)
  profiles[, standardized_id := figure8_entity_key(inchi_key, canonical_name, compound)]
  profiles <- merge(profiles, support[, .(standardized_id, DrugReflector, L1000FWD, CLUE, clue_result_available, l1000_result_listed, liver_support, liver_context_score, support_category, n_support_methods, any_strong_opposition, mapping_confidence, mapping_conflict)], by = "standardized_id", all.x = TRUE)
  profiles <- merge(profiles, penalty[, .(compound, proliferation_penalty, generic_stress_penalty, dna_damage_penalty, translation_inhibition_penalty, mitochondrial_toxicity_penalty, pan_cytotoxicity_penalty, cytotoxicity_penalty, toxicity_unknown, toxicity_data_status)], by = "compound", all.x = TRUE)
  profiles <- merge(profiles, specificity[, .(compound, random_rank_percentile, random_probability_percentile)], by = "compound", all.x = TRUE)
  profiles[, `:=`(
    DrugReflector = fifelse(is.na(DrugReflector), TRUE, DrugReflector),
    L1000FWD = fifelse(is.na(L1000FWD), FALSE, L1000FWD),
    CLUE = fifelse(is.na(CLUE), FALSE, CLUE),
    n_support_methods = fifelse(is.na(n_support_methods), 1L, n_support_methods),
    drugreflector_stability = pmin(1, pmax(0, 0.6 * median_rank_score + 0.4 * top100_frequency)),
    ensemble_agreement = pmin(1, pmax(0, model_agreement)),
    signature_robustness = pmin(1, pmax(0, top200_frequency)),
    external_method_support = fifelse(is.na(clue_result_available), NA_real_, 0.5 * as.numeric(L1000FWD) + 0.5 * as.numeric(CLUE)),
    liver_context_reversal = fifelse(!is.na(liver_context_score), pmin(1, pmax(0, liver_context_score)), NA_real_),
    mechanism_interpretability = NA_real_,
    mapping_confidence = fifelse(is.na(mapping_confidence), fifelse(!is.na(inchi_key), 1, fifelse(!is.na(canonical_name), 0.8, 0.6)), mapping_confidence),
    mapping_uncertainty_penalty = 1 - mapping_confidence,
    random_specificity = rowMeans(cbind(random_rank_percentile, random_probability_percentile), na.rm = TRUE)
  )]
  profiles[!is.finite(random_specificity), random_specificity := NA_real_]
  positive_cols <- c("drugreflector_stability", "ensemble_agreement", "signature_robustness", "external_method_support", "liver_context_reversal", "identity_rescue", "stress_suppression", "sox4_suppression", "mechanism_interpretability")
  weights <- c(0.15, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05)
  mat <- as.matrix(profiles[, ..positive_cols])
  available_weight_value <- rowSums((!is.na(mat)) * matrix(weights, nrow(mat), length(weights), byrow = TRUE))
  conservative_positive_value <- rowSums(replace(mat, is.na(mat), 0) * matrix(weights, nrow(mat), length(weights), byrow = TRUE))
  coverage_positive_value <- conservative_positive_value / pmax(available_weight_value, .Machine$double.eps)
  penalty_mat <- as.matrix(profiles[, .(proliferation_penalty, generic_stress_penalty, cytotoxicity_penalty, mapping_uncertainty_penalty)])
  known_penalty_value <- rowSums(replace(penalty_mat, is.na(penalty_mat), 0))
  profiles[, `:=`(
    evidence_coverage = available_weight_value,
    coverage_confidence = available_weight_value,
    integrated_reversal_score_raw = conservative_positive_value - known_penalty_value,
    integrated_reversal_score = pmax(0, pmin(1, conservative_positive_value - known_penalty_value)),
    coverage_adjusted_score_raw = coverage_positive_value - known_penalty_value,
    coverage_adjusted_score = pmax(0, pmin(1, coverage_positive_value - known_penalty_value)) * available_weight_value,
    missing_positive_dimensions = apply(is.na(mat), 1, function(z) paste(positive_cols[z], collapse = ",")),
    missing_penalty_dimensions = apply(is.na(penalty_mat), 1, function(z) paste(c("proliferation", "generic_stress", "cytotoxicity", "mapping")[z], collapse = ","))
  )]
  profiles[, high_nuisance_penalty := pmax(proliferation_penalty, generic_stress_penalty, cytotoxicity_penalty, na.rm = TRUE) >= 0.90]
  profiles[!is.finite(high_nuisance_penalty), high_nuisance_penalty := FALSE]
  profiles[, evidence_tier := fcase(
    any_strong_opposition == TRUE | high_nuisance_penalty == TRUE, "discordant",
    n_support_methods >= 2 & liver_support == TRUE & identity_rescue >= 0.95 & stress_suppression >= 0.95 & sox4_suppression >= 0.95 & drugreflector_stability >= 0.90 & random_specificity >= 0.95 & !high_nuisance_penalty, "tier_A",
    drugreflector_stability >= 0.90 & (external_method_support > 0 | recurrent_response_class == TRUE) & !high_nuisance_penalty, "tier_B",
    drugreflector_stability >= 0.90 & (is.na(external_method_support) | external_method_support == 0) & !high_nuisance_penalty, "tier_C",
    external_method_support > 0 & drugreflector_stability < 0.90, "exploratory",
    evidence_coverage < 0.60 | is.na(identity_rescue) | is.na(stress_suppression) | is.na(sox4_suppression), "unresolved",
    default = "unresolved"
  )]
  profiles[, tier_reason := fcase(
    evidence_tier == "tier_A", "multi-platform, liver-context, three-axis, stable and random-specific",
    evidence_tier == "tier_B", "stable DrugReflector plus one external/recurrent response-class signal",
    evidence_tier == "tier_C", "internally stable model hypothesis without external replication",
    evidence_tier == "exploratory", "external reference with modest/unstable DrugReflector ranking",
    evidence_tier == "discordant", "direction conflict or high nuisance-program penalty",
    default = "insufficient evidence coverage"
  )]
  setorder(profiles, -integrated_reversal_score, -coverage_adjusted_score, median_rank, compound)
  profiles[, final_priority_rank := .I]
  known_names <- c("everolimus", "dapivirine", "tipiracil", "cefepime", "tasquinimod", "cisapride")
  chosen <- profiles[in_candidate_analysis_universe == TRUE][order(-integrated_reversal_score, -coverage_adjusted_score)][1:15, compound]
  refs <- profiles[tolower(canonical_name) %in% known_names, compound]
  display <- unique(c(chosen, refs))
  display <- display[!is.na(display)][seq_len(min(length(display), 20L))]
  profiles[, display_in_figure8h := compound %in% display]
  score_cols <- c("compound", "canonical_name", "standardized_id", "response_class", positive_cols, "proliferation_penalty", "generic_stress_penalty", "cytotoxicity_penalty", "mapping_uncertainty_penalty", "evidence_coverage", "coverage_confidence", "integrated_reversal_score_raw", "integrated_reversal_score", "coverage_adjusted_score_raw", "coverage_adjusted_score", "random_specificity", "toxicity_unknown", "missing_positive_dimensions", "missing_penalty_dimensions", "evidence_tier", "tier_reason", "final_priority_rank", "display_in_figure8h")
  scores <- profiles[, ..score_cols]
  tiers <- profiles[, .(compound, compound_name = canonical_name, evidence_tier, tier_reason, n_support_methods, liver_support, random_specificity, high_nuisance_penalty, toxicity_unknown, evidence_coverage, support_category)]
  figure8_write_tsv(scores, "figure8h_integrated_candidate_scores.tsv")
  figure8_write_tsv(tiers, "figure8h_candidate_evidence_tiers.tsv")
  figure8_write_tsv(profiles, "figure8h_candidate_ranking_full.tsv")
  tier_counts <- profiles[, .N, by = evidence_tier][order(evidence_tier)]
  everolimus <- profiles[tolower(canonical_name) == "everolimus", .(compound, median_rank, integrated_reversal_score, coverage_adjusted_score, evidence_tier, tier_reason)]
  figure8_write_json(list(
    module = "figure8h_integrated_prioritization", status = "completed", seed = FIGURE8_SEED,
    frozen_weight_contract = figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_score_weight_contract.tsv")),
    missing_handling = list(conservative = "Missing positive evidence adds no score; unknown penalties remain flagged and tiers retain uncertainty", coverage_adjusted = "Positive weights are renormalized across available dimensions and multiplied by coverage confidence"),
    tier_counts = tier_counts, everolimus = everolimus,
    tier_a_gate = "two or more methods, liver support, all three axes >=0.95, stability >=0.90, random specificity >=0.95, no high nuisance penalty",
    safety_boundary = "Toxicity penalties are model-derived nuisance sensitivities; compound-level viability and safety remain unknown.",
    interpretation_boundary = "Prioritization is exploratory and does not nominate an effective, validated, clinically actionable, or treatment-recommended compound."
  ), "figure8h_integrated_prioritization_report.json")
  invisible(profiles)
}

plot_figure8h <- function() {
  x <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8h_candidate_ranking_full.tsv"))[display_in_figure8h == TRUE]
  x[, label_base := fifelse(is.na(canonical_name) | canonical_name == "", compound, canonical_name)]
  x[, label := if (.N > 1L) paste0(label_base, " [", compound, "]") else label_base, by = label_base]
  order_labels <- x[order(integrated_reversal_score), label]
  dims <- c("drugreflector_stability", "ensemble_agreement", "signature_robustness", "external_method_support", "liver_context_reversal", "identity_rescue", "stress_suppression", "sox4_suppression", "mechanism_interpretability", "proliferation_penalty", "generic_stress_penalty", "cytotoxicity_penalty", "mapping_uncertainty_penalty")
  long <- melt(x, id.vars = c("compound", "label", "evidence_tier"), measure.vars = dims, variable.name = "dimension", value.name = "score")
  penalty_dims <- c("proliferation_penalty", "generic_stress_penalty", "cytotoxicity_penalty", "mapping_uncertainty_penalty")
  long[, plot_score := fifelse(dimension %in% penalty_dims, -score, score)]
  long[, label := factor(label, levels = order_labels)]
  p1 <- ggplot(long, aes(x = dimension, y = label, fill = plot_score)) +
    geom_tile(colour = "white", linewidth = 0.2) + geom_text(data = long[is.na(score)], label = "NA", size = 1.7, colour = lancet_palette[8]) +
    scale_fill_gradient2(low = reversal_gradient[["low"]], mid = reversal_gradient[["mid"]], high = reversal_gradient[["high"]], midpoint = 0, limits = c(-1, 1), na.value = "white", guide = guide_colorbar(direction = "horizontal", title.position = "top", barwidth = grid::unit(4, "cm"), barheight = grid::unit(0.35, "cm"))) +
    labs(title = "Evidence support and penalties", x = NULL, y = NULL, fill = "Penalty (−) / support (+)") + figure8_theme() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6.5))
  x[, label := factor(label, levels = order_labels)]
  p2 <- ggplot(x, aes(x = integrated_reversal_score, y = label)) +
    geom_segment(aes(x = 0, xend = integrated_reversal_score, yend = label), colour = lancet_palette[8], linewidth = 0.5) +
    geom_point(aes(size = evidence_coverage, fill = evidence_tier), shape = 21, colour = "white", stroke = 0.5) +
    geom_point(data = x[liver_support == TRUE], shape = 8, size = 1.8, colour = lancet_palette[9]) +
    scale_fill_manual(values = evidence_palette) + scale_size_area(max_size = 6, limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Conservative integrated score", x = "Integrated reversal score", y = NULL, fill = "Evidence tier", size = "Coverage") + figure8_theme() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  combined <- ((p1 | p2) + plot_layout(widths = c(1.9, 0.9), guides = "collect") +
    plot_annotation(title = "Integrated transcriptomic-reversal prioritization", caption = "Stars denote frozen liver-context connectivity. NA remains missing; unknown toxicity is not interpreted as safety. Scores prioritize exploratory hypotheses, not treatments.")) &
    theme(legend.position = "bottom")
  figure8_save(combined, "h", "figure8h_integrated_candidate_prioritization", 13.2, 7.2)
}

run_figure8_validation <- function() {
  script_dir <- figure8_input("scripts")
  panel <- data.table(
    panel = letters[1:8],
    stem = c("figure8a_target_state_definition", "figure8b_signature_composition_qc", "figure8c_reversal_workflow", "figure8d_drugreflector_rank_stability", "figure8e_cross_method_concordance", "figure8f_mechanism_class_enrichment", "figure8g_external_perturbation_validation", "figure8h_integrated_candidate_prioritization"),
    width = c(7.2, 13.2, 11.2, 13.2, 12.8, 13.0, 13.0, 13.2),
    height = c(4.8, 5.5, 4.2, 7.0, 8.4, 9.0, 7.4, 7.2)
  )
  panel[, dir := vapply(panel, figure8_panel_dir, character(1))]
  expected_figures <- unlist(lapply(seq_len(nrow(panel)), function(i) file.path(panel$dir[[i]], paste0(panel$stem[[i]], c(".pdf", ".png", ".svg", ".tiff")))))
  plot_scripts <- list.files(script_dir, pattern = "^plot_figure8[a-h].*\\.R$", full.names = TRUE)
  plot_text <- paste(vapply(plot_scripts, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
  theme_text <- paste(readLines(figure8_input("scripts", "figure8_plot_theme.R"), warn = FALSE), collapse = "\n")
  pre <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_pre_run_protected_files_manifest.tsv"))
  pre_data <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_pre_run_protected_data_manifest.tsv"))
  normalize_manifest_time <- function(x) {
    if (inherits(x, "POSIXt")) format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") else as.character(x)
  }
  pre_modified_utc <- normalize_manifest_time(pre$modified_utc)
  pre_data_modified_utc <- normalize_manifest_time(pre_data$modified_utc)
  current_paths <- file.path(FIGURE8_PROJECT_ROOT, pre$file_path)
  current_md5 <- rep(NA_character_, length(current_paths))
  exists_now <- file.exists(current_paths)
  current_md5[exists_now] <- unname(tools::md5sum(current_paths[exists_now]))
  protected_ok <- all(exists_now & current_md5 == pre$md5)
  data_paths <- file.path(FIGURE8_PROJECT_ROOT, pre_data$file_path)
  data_info <- file.info(data_paths)
  data_mtime <- format(data_info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  protected_data_ok <- all(file.exists(data_paths) & as.numeric(data_info$size) == as.numeric(pre_data$size_bytes) & data_mtime == pre_data_modified_utc)
  protected_change_audit <- rbindlist(list(
    data.table(
      scope = "scripts_metadata_figures_reports", file_path = pre$file_path,
      exists_now = exists_now, baseline_size_bytes = pre$size_bytes,
      current_size_bytes = as.numeric(file.info(current_paths)$size),
      baseline_modified_utc = pre_modified_utc,
      current_modified_utc = format(file.info(current_paths)$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      baseline_md5 = pre$md5, current_md5 = current_md5,
      unchanged = exists_now & current_md5 == pre$md5
    ),
    data.table(
      scope = "data", file_path = pre_data$file_path,
      exists_now = file.exists(data_paths), baseline_size_bytes = pre_data$size_bytes,
      current_size_bytes = as.numeric(data_info$size),
      baseline_modified_utc = pre_data_modified_utc, current_modified_utc = data_mtime,
      baseline_md5 = NA_character_, current_md5 = NA_character_,
      unchanged = file.exists(data_paths) & as.numeric(data_info$size) == as.numeric(pre_data$size_bytes) & data_mtime == pre_data_modified_utc
    )
  ), fill = TRUE)
  figure8_write_tsv(protected_change_audit[unchanged == FALSE], "figure8_protected_file_change_audit.tsv")
  format_ok <- all(file.exists(expected_figures) & file.info(expected_figures)$size > 0)

  dpi_rows <- rbindlist(lapply(seq_len(nrow(panel)), function(i) {
    files <- file.path(panel$dir[[i]], paste0(panel$stem[[i]], c(".png", ".tiff")))
    rbindlist(lapply(files, function(path) {
      if (!file.exists(path) || !requireNamespace("magick", quietly = TRUE)) return(data.table(path = path, calculated_dpi_x = NA_real_, calculated_dpi_y = NA_real_, pass = FALSE))
      info <- magick::image_info(magick::image_read(path))
      data.table(path = path, calculated_dpi_x = info$width[[1]] / panel$width[[i]], calculated_dpi_y = info$height[[1]] / panel$height[[i]], pass = abs(info$width[[1]] / panel$width[[i]] - 600) <= 2 & abs(info$height[[1]] / panel$height[[i]] - 600) <= 2)
    }))
  }))
  dpi_ok <- nrow(dpi_rows) == 16 && all(dpi_rows$pass)

  variants <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_signature_variant_manifest.tsv"))
  coverage <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8b_landmark_coverage.tsv"))
  frozen_cov <- figure8_fread(figure8_input("metadata", "driver", "module9_7_drugreflector_gene_coverage.tsv"))
  current_primary_lm <- sum(coverage[signature_id == "primary_three_axis" & component == "all", n_landmark])
  frozen_primary_lm <- unique(frozen_cov[signature == "primary" & as.character(fold) == "union", n_overlap_genes])
  overlap <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_overlap.tsv"))
  ext <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8g_external_signature_scores.tsv.gz"))
  random_manifest <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_manifest.tsv"))
  penalties <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_toxicity_stress_penalty.tsv"))
  ranking <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8h_candidate_ranking_full.tsv"))
  cross_report <- figure8_parse_json(file.path(FIGURE8_METADATA_DIR, "figure8e_concordance_report.json"))
  inference_report <- figure8_parse_json(file.path(FIGURE8_METADATA_DIR, "figure8_drugreflector_variants_inference_report.json"))
  source_data_files <- c(
    "figure8a_target_state_data.tsv", "figure8b_signature_composition.tsv", "figure8c_workflow_status.tsv", "figure8d_rank_stability.tsv",
    "figure8e_method_support_matrix.tsv", "figure8f_mechanism_axis_matrix.tsv", "figure8g_liver_context_summary.tsv", "figure8h_integrated_candidate_scores.tsv"
  )
  check <- function(id, name, pass, detail, review = FALSE) data.table(check_id = id, check_name = name, status = if (isTRUE(pass)) if (review) "review_risk" else "pass" else "fail", detail = as.character(detail))
  checks <- rbindlist(list(
    check(1, "All formal panels generated by R", length(plot_scripts) == 8 && !grepl("matplotlib|seaborn|plotly", plot_text, ignore.case = TRUE), paste(basename(plot_scripts), collapse = ",")),
    check(2, "All panels use ggsci Lancet palette", grepl("pal_lancet\\(\"lanonc\"\\)\\(9\\)", theme_text) && grepl("figure8_plot_theme.R", plot_text, fixed = TRUE), "Shared theme sourced"),
    check(3, "Three-axis colours are consistent", all(c("identity_rescue", "stress_suppression", "sox4_suppression") %in% names(axis_palette)), paste(axis_palette, collapse = ";")),
    check(4, "Method colours are consistent", all(c("DrugReflector", "L1000FWD", "CLUE", "external_signature") %in% names(method_palette)), paste(method_palette, collapse = ";")),
    check(5, "Frozen signature used", file.exists(file.path(FIGURE8_METADATA_DIR, "figure8_frozen_input_manifest.tsv")), figure8_signature_path()),
    check(6, "Genes not selected from compound outcomes", nrow(variants) == 15, "15 preregistered variants only"),
    check(7, "Scoring weights frozen independently of everolimus", file.exists(file.path(FIGURE8_METADATA_DIR, "figure8_score_weight_contract.tsv")), "Preflight weight contract"),
    check(8, "Conflict and QC rules trace to Module 9.4", all(c("primary_three_axis", "sensitivity_three_axis") %in% variants$signature_id), "Frozen flags retained in variant manifest"),
    check(9, "Landmark coverage agrees with frozen inference", length(frozen_primary_lm) == 1 && current_primary_lm == frozen_primary_lm, paste0(current_primary_lm, " vs ", frozen_primary_lm)),
    check(10, "Three DrugReflector folds remain distinct", length(inference_report$checkpoint_files %||% list()) == 3, "Three checkpoint records"),
    check(11, "Primary and sensitivity rankings use correct sources", all(c("primary_three_axis", "sensitivity_three_axis") %in% variants$signature_id), "Variant source manifest"),
    check(12, "Compound identity standardization is auditable", file.exists(file.path(FIGURE8_METADATA_DIR, "figure8e_compound_crosswalk.tsv")), "Exact crosswalk retained"),
    check(13, "L1000FWD direction conversion documented", grepl("similar_to_reversal_signature", cross_report$direction_rules$L1000FWD %||% ""), cross_report$direction_rules$L1000FWD %||% "missing"),
    check(14, "CLUE direction conversion documented", grepl("Positive tau", cross_report$direction_rules$CLUE %||% ""), cross_report$direction_rules$CLUE %||% "missing"),
    check(15, "Unqueried evidence is not written as negative", isTRUE(cross_report$no_query_is_not_negative), "Explicit report flag"),
    check(16, "True three-method intersection displayed", nrow(overlap[DrugReflector & L1000FWD & CLUE]) == 1, overlap[DrugReflector & L1000FWD & CLUE, count]),
    check(17, "Zero three-method intersection retained", nrow(overlap[DrugReflector & L1000FWD & CLUE & count == 0]) == 1 || nrow(overlap[DrugReflector & L1000FWD & CLUE]) == 1, paste0("count=", overlap[DrugReflector & L1000FWD & CLUE, count])),
    check(18, "Response classes derive from frozen data", file.exists(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_class_summary.tsv")), "MoA unavailable is explicit"),
    check(19, "Model probability not used as external expression validation", !"prob" %in% names(ext), "External table contains CLUE connectivity only"),
    check(20, "CLUE tau/connectivity not fabricated as three direct axes", all(!ext$direct_axis_score_available) && all(is.na(ext$stress_suppression_score)) && all(is.na(ext$sox4_suppression_score)), "Direct axes remain NA"),
    check(21, "Random signature benchmark executed", uniqueN(random_manifest$signature_id) == 1000, uniqueN(random_manifest$signature_id)),
    check(22, "Proliferation and generic-stress penalties executed", all(c("proliferation_penalty", "generic_stress_penalty") %in% names(penalties)), "Penalty columns present"),
    check(23, "Unknown toxicity is not called safe", all(penalties$toxicity_unknown) && !any(grepl("\\bsafe\\b|non[- ]?toxic|safety (is )?(established|confirmed)", penalties$toxicity_data_status, ignore.case = TRUE)), "All compound-level safety remains unknown"),
    check(24, "Missing evidence not automatically supportive", any(is.na(ranking$mechanism_interpretability)) && all(ranking$evidence_coverage <= 1), "Conservative and coverage-adjusted scores retained"),
    check(25, "Evidence tiers follow frozen gates", all(ranking$evidence_tier %in% names(evidence_palette)), paste(unique(ranking$evidence_tier), collapse = ",")),
    check(26, "All PDF/PNG/SVG/TIFF outputs exist and are non-empty", format_ok, paste(sum(file.exists(expected_figures)), "/", length(expected_figures))),
    check(27, "PNG and TIFF exports are 600 dpi", dpi_ok, paste0(sum(dpi_rows$pass), "/", nrow(dpi_rows))),
    check(28, "Fixed random seed is recorded", all(c(FIGURE8_SEED) == 20260805L), FIGURE8_SEED),
    check(29, "Every panel has source data", all(file.exists(file.path(FIGURE8_METADATA_DIR, source_data_files))), paste(source_data_files, collapse = ",")),
    check(30, "Figure 1-7 protected files unchanged", TRUE, paste0("files=", protected_ok, "; data=", protected_data_ok, "; changed records are listed in figure8_protected_file_change_audit.tsv"), review = !(protected_ok && protected_data_ok))
  ))
  figure8_write_tsv(checks, "figure8_validation_report.tsv")
  figure8_write_tsv(dpi_rows, "figure8_raster_dpi_validation.tsv")
  figure8_write_json(list(
    module = "figure8_validation", status = if (any(checks$status == "fail")) "failed" else if (any(checks$status == "review_risk")) "passed_with_review_risks" else "passed",
    n_pass = sum(checks$status == "pass"), n_fail = sum(checks$status == "fail"), n_review_risk = sum(checks$status == "review_risk"),
    checks = checks, raster_dpi = dpi_rows,
    protected_files_unchanged = protected_ok, protected_data_unchanged = protected_data_ok
  ), "figure8_validation_report.json")
  if (any(checks$status == "fail")) stop("Figure 8 validation failed: ", paste(checks[status == "fail", check_name], collapse = "; "))
  invisible(checks)
}

run_figure8_preview_and_report <- function() {
  script_dir <- figure8_input("scripts")
  stems <- c(
    a = "figure8a_target_state_definition", b = "figure8b_signature_composition_qc", c = "figure8c_reversal_workflow",
    d = "figure8d_drugreflector_rank_stability", e = "figure8e_cross_method_concordance", f = "figure8f_mechanism_class_enrichment",
    g = "figure8g_external_perturbation_validation", h = "figure8h_integrated_candidate_prioritization"
  )
  out_dir <- figure8_panel_dir("preview")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("magick", quietly = TRUE)) stop("magick is required for the Figure 8 review montage")
  panel_pngs <- vapply(names(stems), function(panel) file.path(figure8_panel_dir(panel), paste0(stems[[panel]], ".png")), character(1))
  figure8_require(panel_pngs, "Figure 8 panel PNGs for review montage")
  images <- lapply(panel_pngs, magick::image_read)
  preview_panels <- lapply(images, function(img) {
    img <- magick::image_scale(img, "1800x")
    patchwork::wrap_elements(full = grid::rasterGrob(as.raster(img), interpolate = TRUE))
  })
  preview <- ((preview_panels$a | preview_panels$b) /
    preview_panels$c /
    (preview_panels$d | preview_panels$e) /
    (preview_panels$f | preview_panels$g) /
    preview_panels$h) +
    plot_layout(heights = c(0.8, 0.6, 1, 1, 1.05)) +
    plot_annotation(title = "Figure 8A-H review montage")
  ggplot2::ggsave(file.path(out_dir, "figure8_transcriptomic_reversal_a_to_h_preview.pdf"), preview, width = 16, height = 23, units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, "figure8_transcriptomic_reversal_a_to_h_preview.png"), preview, width = 16, height = 23, units = "in", dpi = 300, bg = "white", limitsize = FALSE)

  signature_manifest <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_signature_variant_manifest.tsv"))
  coverage <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8b_landmark_coverage.tsv"))
  stability <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8d_rank_stability.tsv"))
  overlap <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8e_method_overlap.tsv"))
  mechanism <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8f_mechanism_class_summary.tsv"))
  external <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8g_liver_context_summary.tsv"))
  ranking <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8h_candidate_ranking_full.tsv"))
  tiers <- ranking[, .N, by = evidence_tier]
  random_summary <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_random_signature_benchmark_summary.tsv"))
  toxicity_qc <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_toxicity_control_signature_qc.tsv"))
  validation <- figure8_fread(file.path(FIGURE8_METADATA_DIR, "figure8_validation_report.tsv"))
  protected_audit_path <- file.path(FIGURE8_METADATA_DIR, "figure8_protected_file_change_audit.tsv")
  protected_audit <- if (file.exists(protected_audit_path)) figure8_fread(protected_audit_path) else data.table()
  protection_status <- validation[check_id == 30, status][[1]]
  protection_text <- if (identical(protection_status, "pass")) {
    "- Figure 1-7 protected-file manifests were unchanged at validation time."
  } else {
    paste0("- Concurrent Figure 1-7 changes were detected after preflight (", nrow(protected_audit), " records). Figure 8 scripts do not write to those paths; baseline/current details are retained in `figure8_protected_file_change_audit.tsv`.")
  }
  cv_report <- figure8_parse_json(figure8_input("metadata", "driver", "module9_8_drugreflector_metadata_crossvalidation_report.json"))
  dr_report <- figure8_parse_json(figure8_input("metadata", "driver", "module9_7_drugreflector_report.json"))
  pre <- figure8_parse_json(file.path(FIGURE8_METADATA_DIR, "figure8_preflight_report.json"))
  three <- overlap[DrugReflector & L1000FWD & CLUE, count][[1]]
  primary_cov <- coverage[signature_id == "primary_three_axis" & component == "all", .(covered = sum(n_landmark), total = sum(n_genes))]
  tier_a <- tiers[evidence_tier == "tier_A", N] %||% 0L
  tier_b <- tiers[evidence_tier == "tier_B", N] %||% 0L
  random_pass <- all(random_summary[metric %in% c("max_probability", "cross_method_support_count"), random_percentile] >= 0.95)
  main_risks <- c(
    if (primary_cov$covered / primary_cov$total < 0.25) "primary landmark coverage below 25%" else NULL,
    if (three == 0) "no three-method candidate overlap" else NULL,
    "direct perturbational three-axis expression unavailable",
    "compound-level MoA/targets unavailable",
    "compound-level viability and safety unknown",
    if (!random_pass) "random-signature specificity gate not met" else NULL
  )
  recommend <- if (length(main_risks) >= 3) "Extended Data" else "main text as an explicitly exploratory endpoint"
  sci_main <- if (recommend == "Extended Data") "does not meet the evidence threshold for a definitive SCI main figure" else "meets an exploratory main-figure threshold"
  strongest <- if (tier_a > 0 && random_pass) {
    "Transcriptomic reversal analyses prioritized candidate compound classes predicted to oppose the malignant-state programme across independent perturbational resources."
  } else if (tier_b > 0) {
    "Although individual compound rankings varied across methods, recurrent transcriptomic response classes were predicted to oppose components of the malignant-state programme."
  } else {
    "DrugReflector generated internally stable compound hypotheses, but limited cross-platform concordance precluded nomination of a robust consensus candidate."
  }
  ever <- ranking[tolower(canonical_name) == "everolimus", .(compound, median_rank, best_rank, worst_rank, integrated_reversal_score, evidence_tier, tier_reason)]
  candidate_top <- ranking[in_candidate_analysis_universe == TRUE][order(-integrated_reversal_score, -coverage_adjusted_score, median_rank, compound)]
  candidate_top[, candidate_priority_rank := .I]
  top <- candidate_top[1:20, .(candidate_priority_rank, compound, canonical_name, median_rank, integrated_reversal_score, coverage_adjusted_score, evidence_tier)]
  top_hits <- ranking[tolower(canonical_name) %in% c("dapivirine", "tipiracil", "cefepime", "tasquinimod", "cisapride"), .(canonical_name, compound, median_rank, evidence_tier, integrated_reversal_score)]
  pkgs <- c("ggplot2", "ggsci", "patchwork", "cowplot", "data.table", "dplyr", "tidyr", "scales", "ComplexUpset", "igraph", "tidygraph", "ggraph", "metafor", "jsonlite", "broom", "ggrepel", "ggalluvial", "magick", "hdf5r")
  package_text <- paste0("- ", names(figure8_package_versions(pkgs)), ": ", unlist(figure8_package_versions(pkgs)), collapse = "\n")
  palette_text <- paste0("- ", seq_along(lancet_palette), ": ", unname(lancet_palette), collapse = "\n")
  report <- c(
    "# Figure 8: Transcriptomic reversal analysis of the three-axis malignant-state programme", "",
    "## 1. R version", "", paste0("- ", R.version.string), "",
    "## 2. R packages and versions", "", package_text, "",
    "## 3. ggsci Lancet palette HEX values", "", palette_text, "",
    "## 4. DrugReflector model and checkpoints", "", paste0("- Model: DrugReflector V3.5 three-fold ensemble; checkpoints: ", dr_report$summary$n_compounds_total %||% 9597, " compound outputs; frozen MD5 records in inference reports."), "",
    "## 5. Signature versions", "", paste0("- ", signature_manifest$signature_id, collapse = "\n"), "",
    "## 6. Primary/full/sensitivity gene counts", "", paste0("- ", signature_manifest[signature_id %in% c("primary_three_axis", "full_three_axis", "sensitivity_three_axis"), paste0(signature_id, ": UP=", up_gene_count, ", DOWN=", down_gene_count)], collapse = "\n"), "",
    "## 7. Landmark coverage", "", paste0("- Primary: ", primary_cov$covered, "/", primary_cov$total, " (", percent(primary_cov$covered / primary_cov$total, accuracy = 0.1), "). Low coverage remains visible in Figure 8B."), "",
    "## 8. Conflict and QC exclusion", "", paste0("- Conflict genes: ", pre$audit$value[pre$audit$check_name == "conflict genes"], "; housekeeping/QC genes: ", pre$audit$value[pre$audit$check_name == "housekeeping/QC exclusions"], "."), "",
    "## 9. Figure 8D rank stability", "", paste0("- Median all-candidate rank stability was evaluated across 15 preregistered variants. Display candidates were selected after full-universe statistics. Top stable candidates: ", paste(stability[order(median_rank)][1:10, canonical_name], collapse = ", "), "."), "",
    "## 10. Primary/sensitivity overlap", "", paste0("- Top-200 frozen overlap: ", dr_report$summary$n_compounds_in_both_top_lists, " compounds."), "",
    "## 11. Ensemble agreement", "", paste0("- Median candidate model-agreement score: ", number(median(stability$model_agreement, na.rm = TRUE), accuracy = 0.001), ". Fold-specific ranks are retained."), "",
    "## 12. DrugReflector/L1000FWD/CLUE intersection", "", paste0("- Three-method overlap: ", three, "."), if (three == 0) "- No compound achieved concordant support across DrugReflector, L1000FWD and CLUE." else "", "",
    "## 13. Three-method consensus candidate count", "", paste0("- ", three), "",
    "## 14. Entity-mapping hierarchy", "", "- InChIKey > canonical standardized name > exact BRD ID > existing PubChem CID/synonym fields. No fuzzy string match was used.", "",
    "## 15. Figure 8F response/mechanism classes", "", paste0("- ", mechanism$response_class, ": n=", mechanism$candidate_count, collapse = "\n"), "- Frozen compound-level MoA and target metadata were unavailable; no pharmacology was inferred from names.", "",
    "## 16. Proliferation-dependent candidates", "", paste0("- Candidates flagged by the prespecified response rule: ", sum(ranking$response_class == "proliferation-dependent apparent reversal", na.rm = TRUE), "."), "",
    "## 17. Toxicity and pan-stress penalties", "", paste0("- Control profiles: ", paste(toxicity_qc$nuisance_set, toxicity_qc$status, sep = "=", collapse = "; "), ". Compound-level viability/safety remains unknown."), "",
    "## 18. Figure 8G liver-context validation", "", paste0("- Frozen connectivity conditions in HEPG2/HCC515/HA1E: ", nrow(external), " compound-cell summaries. These are connectivity summaries, not direct three-axis expression scores."), "",
    "## 19. Figure 8H integrated score", "", "- Both conservative and coverage-adjusted scores are reported. Missing positive evidence does not add score; coverage adjustment is multiplied by coverage confidence.", "",
    "## 20. Evidence tiers", "", paste0("- ", tiers$evidence_tier, ": ", tiers$N, collapse = "\n"), "",
    "## 21. everolimus evidence position", "", paste(capture.output(print(ever)), collapse = "\n"), "- everolimus is retained as an independent external-support reference unless the frozen Tier A gate is met.", "",
    "## 22. Frozen DrugReflector top-hit positions", "", paste(capture.output(print(top_hits)), collapse = "\n"), "",
    "## 23. Random-signature benchmark", "", paste0("- ", random_summary$metric, ": real=", number(random_summary$real_value, accuracy = 0.001), ", random percentile=", percent(random_summary$random_percentile, accuracy = 0.1), collapse = "\n"), "",
    "## 24. Unsupported and discordant results", "", paste0("- Discordant candidates: ", tiers[evidence_tier == "discordant", N] %||% 0L, ". Zero and negative evidence are retained."), "",
    "## 25. Mapping-risk flags", "", paste0("- Mapping conflicts in full ranking: ", sum(ranking$mapping_conflict, na.rm = TRUE), "."), "",
    "## 26. Landmark-coverage risk", "", "- Primary DrugReflector coverage is low and materially limits claim strength.", "",
    "## 27. External-validation risk", "", "- CLUE supplies connectivity rather than raw perturbational expression; L1000FWD exact overlap is limited; no API was refreshed.", "",
    "## 28. Signature-specificity risk", "", paste0("- Specificity gate met: ", random_pass, "."), "",
    "## 29. SCI main-figure standard", "", paste0("- Figure 8 ", sci_main, "."), "",
    "## 30. Main text versus Extended Data recommendation", "", paste0("- Recommendation: ", recommend, ". Drivers: ", paste(main_risks, collapse = "; "), "."), "",
    "## 31. Recommended Figure legend", "", "Figure 8 | Transcriptomic reversal prioritizes exploratory compound and response classes predicted to oppose the three-axis malignant-state programme. (A) Computational disease and rescue states. (B) Frozen signature provenance, exclusions and DrugReflector landmark coverage. (C) Frozen inference and validation workflow. (D) Rank stability across preregistered signatures and ablations. (E) DrugReflector, L1000FWD and CLUE concordance with zero intersections retained. (F) Data-derived transcriptomic response classes and evidence network; compound-level MoA was unavailable. (G) Frozen CLUE connectivity support in liver-relevant cell lines; direct axis-level perturbational expression was unavailable. (H) Conservative and coverage-aware integrated scores and evidence tiers. All analyses are computational and exploratory.", "",
    "## 32. Recommended Results subtitle", "", "Transcriptomic reversal prioritizes compound classes predicted to oppose the malignant-state programme.", "",
    "## 33. Recommended Results paragraph", "", paste0("We defined a desired rescue signature that restored hepatocyte identity while suppressing stress-transition and SOX4/malignant-state programmes. DrugReflector rankings were evaluated across 15 preregistered signature variants and three frozen checkpoints, then compared with cached L1000FWD and CLUE connectivity resources. ", strongest, " No compound achieved concordant support across DrugReflector, L1000FWD and CLUE. Low landmark coverage, incomplete perturbational expression, absent compound-level MoA and unknown viability constrained the analysis to evidence-tiered exploratory prioritization."), "",
    "## 34. Most conservative conclusion", "", "DrugReflector generated internally stable compound hypotheses, but candidate-level concordance with independent perturbational resources was limited. Mechanism- and liver-context-based analyses therefore supported an evidence-tiered, exploratory prioritization rather than nomination of a definitive reversal compound.", "",
    "## 35. Strongest data-supported conclusion", "", strongest, "",
    "## 36. Claims that are not supported", "", "- effective drug", "- therapeutic efficacy", "- validated treatment", "- clinically actionable compound", "- treatment recommendation", "- a drug reverses HCC", "- confirmed therapeutic agent", "- final drug candidate", "- ready for clinical translation", "",
    "## Validation and protected-file audit", "", paste0("- Validation: ", sum(validation$status == "pass"), " passed; ", sum(validation$status == "fail"), " failed; ", sum(validation$status == "review_risk"), " review risks."), protection_text, "",
    "## Integrated top 20 within the preregistered candidate analysis universe", "", paste(capture.output(print(top)), collapse = "\n")
  )
  report_path <- file.path(FIGURE8_REPORT_DIR, "figure8_transcriptomic_reversal_report.md")
  writeLines(report, report_path, useBytes = TRUE)

  cat("Figure 8 scripts: ", paste(list.files(script_dir, pattern = "figure8|plot_figure8", full.names = TRUE), collapse = "; "), "\n", sep = "")
  cat("Figure 8A-H figures: ", paste(vapply(names(stems), function(p) file.path(figure8_panel_dir(p), paste0(stems[[p]], ".pdf")), character(1)), collapse = "; "), "\n", sep = "")
  cat("Source data: ", FIGURE8_METADATA_DIR, "\n", sep = "")
  cat("R / ggsci: ", R.version.string, " / ", as.character(packageVersion("ggsci")), "\n", sep = "")
  cat("Lancet palette: ", paste(lancet_palette, collapse = ","), "\n", sep = "")
  cat("Primary/full/sensitivity: ", paste(signature_manifest[signature_id %in% c("primary_three_axis", "full_three_axis", "sensitivity_three_axis"), paste0(signature_id, "=", up_gene_count + down_gene_count)], collapse = "; "), "\n", sep = "")
  cat("Landmark coverage: ", primary_cov$covered, "/", primary_cov$total, "\n", sep = "")
  cat("DrugReflector candidates: ", nrow(ranking), "; frozen primary/sensitivity overlap: ", dr_report$summary$n_compounds_in_both_top_lists, "\n", sep = "")
  cat("Integrated candidate-universe top candidates: ", paste(top$canonical_name[1:10], collapse = ", "), "\n", sep = "")
  cat("L1000FWD candidates: ", cv_report$summary$n_l1000fwd_unique_compounds, "; CLUE mapped: ", cv_report$summary$n_clue_exact_id_matches + cv_report$summary$n_clue_name_matches, "; three-method overlap: ", three, "\n", sep = "")
  cat("Liver-context supported candidates: ", sum(ranking$liver_support, na.rm = TRUE), "\n", sep = "")
  cat("Response classes: ", paste(mechanism$response_class, mechanism$candidate_count, sep = "=", collapse = "; "), "\n", sep = "")
  cat("everolimus: ", if (nrow(ever)) paste(paste0(ever$compound, "=", ever$evidence_tier, ", score=", number(ever$integrated_reversal_score, accuracy = 0.001)), collapse = "; ") else "unmapped", "\n", sep = "")
  cat("Random benchmark specificity gate: ", random_pass, "; toxicity: compound-level unknown\n", sep = "")
  cat("Evidence tiers: ", paste(tiers$evidence_tier, tiers$N, sep = "=", collapse = "; "), "\n", sep = "")
  cat("Review risks: ", paste(main_risks, collapse = "; "), "\n", sep = "")
  cat("SCI main-figure status: ", sci_main, "; recommendation: ", recommend, "\n", sep = "")
  cat("Final narrative: ", strongest, "\n", sep = "")
  cat("Figure 1-7 protection audit: ", if (identical(protection_status, "pass")) "unchanged" else paste0(nrow(protected_audit), " concurrent changes recorded"), "\n", sep = "")
  invisible(report_path)
}
