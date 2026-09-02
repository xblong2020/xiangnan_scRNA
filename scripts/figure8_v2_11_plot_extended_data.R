#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_theme.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(tidygraph)
  library(ggraph)
  library(ComplexUpset)
})

figure8_v2_plot_extended <- function() {
  evidence <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_integrated_candidate_evidence.tsv"))

  predictions <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_drugreflector_variant_predictions.tsv.gz"))
  display <- evidence[candidate_analysis_universe == TRUE][order(candidate_priority_rank)][1:50, compound]
  heat <- predictions[compound %in% display]
  labels <- evidence[compound %in% display, .(compound, label = make.unique(canonical_name))]
  heat <- merge(heat, labels, by = "compound")
  heat[, rank_signal := -log10(rank_1based / 9597)]
  heat[, label := factor(label, levels = rev(labels$label))]
  extended_1 <- ggplot(heat, aes(x = signature_id, y = label, fill = rank_signal)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_gradientn(colours = c("#F7F7F7", figure8_v2_lancet[[6]], figure8_v2_lancet[[2]])) +
    labs(title = "Extended Data 8-1 | All 16 signature definitions", x = NULL, y = NULL, fill = expression(-log[10](rank/N))) +
    figure8_v2_theme() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.5))
  figure8_v2_save_plot(extended_1, "extended_data_8_1_signature_rank_heatmap", 183, 145)

  cross <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_concordance.tsv"))
  upset_data <- as.data.frame(cross[, .(DrugReflector = v2_drugreflector_support, L1000FWD, CLUE)])
  p_upset <- ComplexUpset::upset(upset_data, intersect = c("DrugReflector", "L1000FWD", "CLUE"), min_size = 0) +
    ggtitle("Cross-framework entity intersections")
  combos <- CJ(DrugReflector = c(FALSE, TRUE), L1000FWD = c(FALSE, TRUE), CLUE = c(FALSE, TRUE))
  combo_counts <- cross[, .N, by = .(DrugReflector = v2_drugreflector_support, L1000FWD, CLUE)]
  combos <- merge(combos, combo_counts, by = c("DrugReflector", "L1000FWD", "CLUE"), all.x = TRUE)
  combos[is.na(N), N := 0L]
  combos[, combination := paste0(ifelse(DrugReflector, "DR", ""), ifelse(L1000FWD, "+L1000", ""), ifelse(CLUE, "+CLUE", ""))]
  combos[combination == "", combination := "none"]
  p_combo <- ggplot(combos, aes(x = N, y = reorder(combination, N), fill = N > 0)) +
    geom_col() + geom_text(aes(label = N), hjust = -0.1, size = 2.2) +
    scale_fill_manual(values = c(`TRUE` = figure8_v2_lancet[[1]], `FALSE` = "white")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "All combinations (zero retained)", x = "Standardized entities", y = NULL) + figure8_v2_theme() + theme(legend.position = "none")
  extended_2 <- p_upset / p_combo + plot_annotation(title = "Extended Data 8-2 | Full DrugReflector/L1000FWD/CLUE overlap")
  figure8_v2_save_plot(extended_2, "extended_data_8_2_cross_framework_upset", 183, 160)

  nodes <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8f_evidence_network_nodes.tsv"))
  edges <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8f_evidence_network_edges.tsv"))
  graph <- tbl_graph(nodes = nodes[, .(name, label, node_type)], edges = edges[, .(from, to, edge_type)], directed = FALSE)
  extended_3 <- ggraph(graph, layout = "stress") +
    geom_edge_link(aes(linetype = edge_type), colour = "#B0B0B0", alpha = 0.75) +
    geom_node_point(aes(fill = node_type), shape = 21, size = 3, colour = "white") +
    geom_node_text(aes(label = label), repel = TRUE, size = 2) +
    scale_fill_manual(values = c(compound = figure8_v2_lancet[[6]], method = figure8_v2_lancet[[1]], response_class = figure8_v2_lancet[[3]])) +
    labs(title = "Extended Data 8-3 | Full transcriptomic response-class network", fill = "Node", edge_linetype = "Edge") +
    theme_void(base_size = 6.5) + theme(legend.position = "bottom", plot.title = element_text(size = 7, face = "bold"))
  figure8_v2_save_plot(extended_3, "extended_data_8_3_response_class_network", 183, 145)

  lincs <- figure8_v2_read_tsv(file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal/figure8g_liver_context_summary.tsv"))
  lincs[, corrected_context := fifelse(cell_line == "HEPG2", "liver-derived / hepatoblastoma-like caveat", fifelse(cell_line == "HCC515", "lung adenocarcinoma / non-liver", fifelse(cell_line == "HA1E", "kidney-derived / non-liver", "other")))]
  lincs_display <- lincs[display_in_main_panel == TRUE]
  extended_4 <- ggplot(lincs_display, aes(x = mean_combined_connectivity, y = reorder(compound_name, mean_combined_connectivity), colour = corrected_context)) +
    geom_vline(xintercept = 0, colour = "#A0A0A0") + geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, linewidth = 0.35) +
    geom_point(size = 2) + facet_wrap(~ corrected_context, scales = "free_y") +
    scale_colour_manual(values = c(`liver-derived / hepatoblastoma-like caveat` = figure8_v2_lancet[[1]], `lung adenocarcinoma / non-liver` = figure8_v2_lancet[[3]], `kidney-derived / non-liver` = figure8_v2_lancet[[4]], other = figure8_v2_lancet[[8]])) +
    labs(title = "Extended Data 8-4 | Frozen LINCS cell-context connectivity with corrected annotations", x = "Mean combined connectivity (descriptive CI)", y = NULL, colour = NULL) + figure8_v2_theme() + theme(legend.position = "bottom")
  figure8_v2_save_plot(extended_4, "extended_data_8_4_lincs_context_connectivity", 183, 145)

  prism_cell <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_prism_primary_candidate_cell_values.tsv.gz"))
  prism_cell[, display_name := factor(brd_core, levels = unique(brd_core))]
  prism_cell[, context_group := fifelse(verified_context == "adult_hepatocellular_carcinoma", "adult HCC", fifelse(grepl("hepatoblastoma", verified_context), "hepatoblastoma-like", "other cancer"))]
  extended_5 <- ggplot(prism_cell, aes(x = factor(depmap_id), y = display_name, fill = primary_sensitivity)) +
    geom_tile() + facet_grid(. ~ context_group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = figure8_v2_lancet[[2]], mid = "#F7F7F7", high = figure8_v2_lancet[[1]], midpoint = 0) +
    labs(title = "Extended Data 8-5 | PRISM primary pan-cancer viability matrix for mapped candidates", x = "Cell models (labels suppressed)", y = "BRD core", fill = "-LFC sensitivity") +
    figure8_v2_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  figure8_v2_save_plot(extended_5, "extended_data_8_5_prism_pan_cancer_matrix", 183, 130)

  null <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_inference_summary.tsv.gz"))
  null_long <- melt(null, id.vars = "signature_id", measure.vars = c("max_probability", "top10_probability_concentration", "top100_probability_concentration", "top_candidate_model_agreement"), variable.name = "metric", value.name = "null_value")
  observed <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_random_specificity_summary.tsv"))
  observed <- observed[metric %in% unique(null_long$metric), .(metric, observed_value, empirical_p_two_sided, specificity_status)]
  extended_6 <- ggplot(null_long, aes(x = null_value)) +
    geom_histogram(bins = 35, fill = figure8_v2_lancet[[8]], colour = "white") +
    geom_vline(data = observed, aes(xintercept = observed_value), colour = figure8_v2_lancet[[2]], linewidth = 0.6) +
    facet_wrap(~ metric, scales = "free", ncol = 2) +
    labs(title = "Extended Data 8-6 | Matched-random null distributions", x = "Null metric", y = "Count") + figure8_v2_theme()
  figure8_v2_save_plot(extended_6, "extended_data_8_6_matched_random_distributions", 183, 130)

  annotation <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_compound_moa_target_annotation.tsv"))
  audit <- annotation[, .(
    `MoA available` = sum(!is.na(curated_MoA)),
    `Target available` = sum(!is.na(curated_targets)),
    `Both available` = sum(!is.na(curated_MoA) & !is.na(curated_targets)),
    `Mapping conflict` = sum(mapping_conflict, na.rm = TRUE),
    `Unavailable` = sum(is.na(curated_MoA) & is.na(curated_targets))
  )]
  audit_long <- melt(audit, variable.name = "audit_class", value.name = "n_compounds")
  extended_7 <- ggplot(audit_long, aes(x = n_compounds, y = reorder(audit_class, n_compounds), fill = audit_class)) +
    geom_col() + geom_text(aes(label = n_compounds), hjust = -0.1, size = 2.3) +
    scale_fill_manual(values = setNames(rep(figure8_v2_lancet[c(1, 3, 4, 2, 8)], length.out = nrow(audit_long)), audit_long$audit_class)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Extended Data 8-7 | Exact entity mapping and curated MoA/target audit", x = "DrugReflector entities", y = NULL) + figure8_v2_theme() + theme(legend.position = "none")
  figure8_v2_save_plot(extended_7, "extended_data_8_7_mapping_moa_audit", 183, 90)

  cat("FIGURE8_V2_EXTENDED_DATA exported=7\n")
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") figure8_v2_plot_extended()

