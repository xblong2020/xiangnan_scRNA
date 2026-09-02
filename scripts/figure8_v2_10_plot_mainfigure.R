#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_theme.R"))
source(file.path("scripts", "figure8_v2_06_moa_network.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggrepel)
  library(tidygraph)
  library(ggraph)
})

figure8_v2_plot_main <- function() {
  figure8_v2_init_dirs()
  evidence <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_integrated_candidate_evidence.tsv"))
  top20 <- evidence[candidate_analysis_universe == TRUE][order(candidate_priority_rank)][1:20]

  state <- data.table(
    programme = factor(c("Identity", "Stress transition", "SOX4/plasticity", "Malignant fate"), levels = rev(c("Identity", "Stress transition", "SOX4/plasticity", "Malignant fate"))),
    malignant = c(-1, 1, 1, 1), rescue = c(1, -1, -1, -1),
    axis = c("axis_A_identity", "axis_B_stress", "axis_C_sox4", "malignant_fate")
  )
  figure8_v2_write_tsv(state, "figure8_v2_panel_a_source.tsv")
  panel_a <- ggplot(state, aes(y = programme)) +
    geom_segment(aes(x = malignant, xend = rescue, yend = programme), colour = "#B8B8B8", linewidth = 0.45) +
    geom_point(aes(x = malignant, fill = axis), shape = 21, size = 3.1, colour = "white", stroke = 0.35) +
    geom_point(aes(x = rescue, fill = axis), shape = 21, size = 3.1, colour = "white", stroke = 0.35) +
    scale_fill_manual(values = figure8_v2_axis_palette) +
    scale_x_continuous(breaks = c(-1, 1), labels = c("Malignant state", "Desired rescue"), limits = c(-1.35, 1.35)) +
    labs(title = "Malignant and desired rescue states", x = NULL, y = NULL,
         caption = "Transcriptomic rescue is not tumour killing or demonstrated phenotypic rescue.") +
    figure8_v2_theme() + theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank())

  score <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_gene_level_rescue_vscore.tsv"))
  balance <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_landmark_axis_balance.tsv"))
  landmark_counts <- data.table(
    category = factor(c("978 model space", "Non-zero scored", "Positive", "Negative", "v1 sparse overlap"), levels = rev(c("978 model space", "Non-zero scored", "Positive", "Negative", "v1 sparse overlap"))),
    count = c(978, sum(abs(score$final_rescue_vscore) > 0), sum(score$final_rescue_vscore > 0), sum(score$final_rescue_vscore < 0), 47),
    denominator = c(978, 978, 978, 978, 300),
    group = c("v2", "v2", "positive", "negative", "v1")
  )
  figure8_v2_write_tsv(landmark_counts, "figure8_v2_panel_b_landmark_source.tsv")
  p_b1 <- ggplot(landmark_counts, aes(x = count / denominator, y = category, fill = group)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = paste0(count, "/", denominator)), hjust = -0.08, size = 1.8) +
    scale_fill_manual(values = c(v2 = figure8_v2_lancet[[1]], positive = figure8_v2_lancet[[1]], negative = figure8_v2_lancet[[2]], zero = figure8_v2_lancet[[8]], v1 = figure8_v2_lancet[[5]])) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1.18)) +
    labs(title = "Continuous landmark-space representation", x = "Fraction of stated universe", y = NULL) +
    figure8_v2_theme() + theme(legend.position = "none")
  p_b2 <- ggplot(balance, aes(x = absolute_mass_fraction, y = axis, fill = axis)) +
    geom_col(width = 0.58) + geom_text(aes(label = scales::percent(absolute_mass_fraction, accuracy = 0.1)), hjust = -0.08, size = 2) +
    scale_fill_manual(values = c(axis_A = figure8_v2_axis_palette[["axis_A_identity"]], axis_B = figure8_v2_axis_palette[["axis_B_stress"]], axis_C = figure8_v2_axis_palette[["axis_C_sox4"]])) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 0.76)) +
    labs(x = "Three-axis absolute score mass", y = NULL, caption = "No severe axis domination by the frozen >70% / effective-axis <1.8 rule.") +
    figure8_v2_theme() + theme(legend.position = "none")
  panel_b <- patchwork::wrap_elements(full = patchwork::patchworkGrob(p_b1 / p_b2 + plot_layout(heights = c(1.25, 0.8))))

  display_c <- top20[1:6, .(canonical_name, compound, median_rank, rank_q1, rank_q3, best_rank, worst_rank, fold_model_agreement)]
  display_c[, label := fifelse(duplicated(canonical_name) | duplicated(canonical_name, fromLast = TRUE), paste0(canonical_name, " [", compound, "]"), canonical_name)]
  display_c[, label := factor(label, levels = rev(label))]
  figure8_v2_write_tsv(display_c, "figure8_v2_panel_c_rank_source.tsv")
  p_c1 <- ggplot(display_c, aes(x = median_rank, y = label)) +
    geom_segment(aes(x = best_rank, xend = worst_rank, yend = label), colour = "#B8B8B8", linewidth = 0.45) +
    geom_segment(aes(x = rank_q1, xend = rank_q3, yend = label), colour = "#333333", linewidth = 1.2) +
    geom_point(aes(fill = fold_model_agreement), shape = 21, size = 3.0, colour = "white", stroke = 0.35) +
    scale_x_log10() + scale_fill_gradient(low = figure8_v2_lancet[[8]], high = figure8_v2_lancet[[1]], limits = c(0, 1)) +
    labs(title = "Internal rank robustness", x = "Rank (median, IQR and range)", y = NULL, fill = "Fold agreement") + figure8_v2_theme() + theme(axis.text.y = element_text(size = 5.3))
  random <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_matched_random_inference_summary.tsv.gz"))
  specificity <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_random_specificity_summary.tsv"))
  observed_top100 <- specificity[metric == "top100_probability_concentration", observed_value]
  observed_p <- specificity[metric == "top100_probability_concentration", empirical_p_two_sided]
  figure8_v2_write_tsv(random[, .(signature_id, top100_probability_concentration)], "figure8_v2_panel_c_null_source.tsv")
  p_c2 <- ggplot(random, aes(x = top100_probability_concentration)) +
    geom_histogram(bins = 35, fill = figure8_v2_lancet[[8]], colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = observed_top100, colour = figure8_v2_lancet[[2]], linewidth = 0.7) +
    annotate("text", x = observed_top100, y = Inf, label = paste0("Observed\nPemp=", format(observed_p, digits = 3)), vjust = 1.2, hjust = -0.08, size = 2.1, colour = figure8_v2_lancet[[2]]) +
    labs(title = "Matched-random specificity", x = "Top-100 probability concentration", y = "2,000 matched nulls") + figure8_v2_theme()
  panel_c <- patchwork::wrap_elements(full = patchwork::patchworkGrob(p_c1 / p_c2 + plot_layout(heights = c(1.25, 0.75))))

  cross <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_cross_framework_concordance.tsv"))
  display_d <- cross[n_support_frameworks >= 2][order(-n_support_frameworks, v2_drugreflector_rank)]
  display_d[, label := fifelse(is.na(canonical_name), standardized_id, canonical_name)]
  display_d[, label := make.unique(label)]
  cross_long <- melt(
    display_d[, .(label, v2_drugreflector_support, L1000FWD, CLUE, l1000_result_available, clue_result_available, external_strong_opposition)],
    id.vars = c("label", "l1000_result_available", "clue_result_available", "external_strong_opposition"),
    measure.vars = c("v2_drugreflector_support", "L1000FWD", "CLUE"), variable.name = "framework", value.name = "support"
  )
  cross_long[, framework := factor(framework, levels = c("v2_drugreflector_support", "L1000FWD", "CLUE"), labels = c("DrugReflector", "L1000FWD", "CLUE"))]
  cross_long[, state := fifelse(support, "support", fifelse(framework == "L1000FWD" & !l1000_result_available | framework == "CLUE" & !clue_result_available, "not available", fifelse(external_strong_opposition, "discordant", "profiled / no support")))]
  cross_long[, label := factor(label, levels = rev(unique(display_d$label)))]
  figure8_v2_write_tsv(cross_long, "figure8_v2_panel_d_source.tsv")
  panel_d <- ggplot(cross_long, aes(x = framework, y = label, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c(support = figure8_v2_lancet[[1]], discordant = figure8_v2_lancet[[2]], `profiled / no support` = "#F2F2F2", `not available` = figure8_v2_lancet[[8]])) +
    labs(title = "Cross-framework transcriptomic corroboration", subtitle = "v1 three-way overlap = 0; v2 continuous-primary overlap = 3", x = NULL, y = NULL, fill = NULL,
         caption = "L1000FWD and CLUE are related LINCS/CMap-derived resources.") + figure8_v2_theme() + theme(legend.position = "bottom")

  target_edges <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_candidate_target_edges.tsv"))
  network <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_network_consistency.tsv"))
  selected_brd <- top20[!is.na(curated_targets), compound][1:min(8, sum(!is.na(top20$curated_targets)))]
  edges <- target_edges[BRD_ID %in% selected_brd]
  candidate_labels <- unique(edges[, .(node = candidate, node_type = "candidate")])
  target_labels <- unique(edges[, .(node = target, node_type = "curated_target")])
  graph_edges <- edges[, .(from = candidate, to = target, edge_type = "curated target")]
  compat <- network[BRD_ID %in% selected_brd & !is.na(compatible_axes)]
  if (nrow(compat)) {
    compat_edges <- compat[, {
      axes <- unlist(strsplit(compatible_axes, ";", fixed = TRUE))
      targets <- figure8_v2_split_targets(curated_targets)
      if (!length(targets) || !length(axes)) NULL else CJ(from = targets, to = axes)[, edge_type := "network-consistent inference"]
    }, by = BRD_ID][, BRD_ID := NULL]
    graph_edges <- rbindlist(list(graph_edges, compat_edges), fill = TRUE)
  }
  axis_nodes <- unique(graph_edges[grepl("^axis_", to), .(node = to, node_type = "axis")])
  graph_nodes <- unique(rbindlist(list(candidate_labels, target_labels, axis_nodes), fill = TRUE), by = "node")
  figure8_v2_write_tsv(graph_edges, "figure8_v2_panel_e_edges.tsv")
  figure8_v2_write_tsv(graph_nodes, "figure8_v2_panel_e_nodes.tsv")
  if (nrow(graph_edges)) {
    graph <- tidygraph::tbl_graph(nodes = graph_nodes, edges = graph_edges, directed = TRUE)
    panel_e <- ggraph(graph, layout = "stress") +
      geom_edge_link(aes(linetype = edge_type), colour = "#A0A0A0", alpha = 0.8, end_cap = circle(2, "mm"), arrow = grid::arrow(length = grid::unit(1.4, "mm"))) +
      geom_node_point(aes(fill = node_type), shape = 21, size = 3.2, colour = "white", stroke = 0.35) +
      geom_node_text(aes(label = node), repel = TRUE, size = 2.1) +
      scale_fill_manual(values = c(candidate = figure8_v2_lancet[[6]], curated_target = figure8_v2_lancet[[1]], axis = figure8_v2_lancet[[3]])) +
      scale_edge_linetype_manual(values = c(`curated target` = "solid", `network-consistent inference` = "dashed")) +
      labs(title = "Curated mechanisms connect candidates to the frozen network", fill = "Node", edge_linetype = "Evidence") +
      theme_void(base_size = 6.5) + theme(legend.position = "none", plot.title = element_text(size = 7, face = "bold"))
  } else panel_e <- ggplot() + annotate("text", x = 0, y = 0, label = "No exact curated target edges for displayed candidates") + labs(title = "Curated mechanisms and targets") + theme_void()

  prism <- figure8_v2_read_tsv(file.path(FIGURE8_V2_METADATA, "figure8_v2_prism_viability.tsv"))
  display_f <- prism[candidate_analysis_universe == TRUE & secondary_dose_response_available == TRUE][order(v2_primary_rank)][1:20]
  enriched <- prism[prism_gate_class == "hcc_liver_enriched"]
  display_f <- unique(rbindlist(list(display_f, enriched), fill = TRUE), by = "BRD_ID")
  display_f[, label := canonical_name]
  figure8_v2_write_tsv(display_f, "figure8_v2_panel_f_source.tsv")
  panel_f <- ggplot(display_f, aes(x = secondary_pan_cancer_median_auc, y = secondary_adult_hcc_median_auc)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#A0A0A0") +
    geom_point(aes(fill = prism_phenotype_class, size = n_secondary_adult_hcc_lines), shape = 21, colour = "white", stroke = 0.4) +
    ggrepel::geom_text_repel(aes(label = label), size = 1.6, max.overlaps = 10, min.segment.length = 0) +
    scale_fill_manual(values = c(hcc_liver_enriched_single_release = figure8_v2_lancet[[3]], no_enriched_support = figure8_v2_lancet[[8]], pan_cancer_activity = figure8_v2_lancet[[2]], unavailable = "white"), na.value = "white") +
    scale_size_area(max_size = 5) +
    labs(title = "Orthogonal PRISM cancer-cell viability", x = "Pan-cancer median AUC", y = "Adult-HCC median AUC", fill = "Phenotype", size = "HCC lines",
         caption = "Lower AUC indicates greater sensitivity; PRISM is not normal-cell safety.") + figure8_v2_theme() + theme(legend.position = "none")

  matrix_cols <- c("DR_score", "robustness_score", "signature_specificity_score", "cross_framework_score", "network_moa_score", "prism_phenotype_score", "nuisance_penalty")
  display_g <- top20[1:12]
  matrix_long <- melt(display_g[, c("canonical_name", "compound", "evidence_tier", "conservative_score", "evidence_coverage", matrix_cols), with = FALSE], id.vars = c("canonical_name", "compound", "evidence_tier", "conservative_score", "evidence_coverage"), measure.vars = matrix_cols, variable.name = "component", value.name = "value")
  matrix_long[, label := fifelse(duplicated(canonical_name) | duplicated(canonical_name, fromLast = TRUE), paste0(canonical_name, " [", compound, "]"), canonical_name)]
  label_order <- unique(matrix_long$label)
  matrix_long[, label := factor(label, levels = rev(label_order))]
  matrix_long[, plot_value := fifelse(component == "nuisance_penalty", -value, value)]
  matrix_long[, component := factor(component, levels = matrix_cols, labels = c("DR", "Robustness", "Specificity", "CMap", "MoA/network", "PRISM", "Nuisance"))]
  figure8_v2_write_tsv(matrix_long, "figure8_v2_panel_g_source.tsv")
  p_g1 <- ggplot(matrix_long, aes(x = component, y = label, fill = plot_value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(data = matrix_long[is.na(value)], label = "NA", size = 1.55, colour = "#999999") +
    scale_fill_gradient2(low = figure8_v2_diverging[["low"]], mid = figure8_v2_diverging[["mid"]], high = figure8_v2_diverging[["high"]], midpoint = 0, limits = c(-1, 1), na.value = "white") +
    labs(title = "Integrated multi-layer evidence prioritizes exploratory hypotheses", x = NULL, y = NULL, fill = "Penalty (-) / support (+)") +
    figure8_v2_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.2), axis.text.y = element_text(size = 5.3), legend.position = "bottom")
  score_points <- unique(matrix_long[, .(label, evidence_tier, conservative_score, evidence_coverage)])
  p_g2 <- ggplot(score_points, aes(x = conservative_score, y = label)) +
    geom_segment(aes(x = 0, xend = conservative_score, yend = label), colour = "#B8B8B8", linewidth = 0.4) +
    geom_point(aes(fill = evidence_tier, size = evidence_coverage), shape = 21, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = figure8_v2_tier_palette) + scale_size_area(max_size = 4.5, limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Conservative evidence score", y = NULL, fill = "Tier", size = "Coverage") + figure8_v2_theme() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  panel_g <- patchwork::wrap_elements(full = patchwork::patchworkGrob((p_g1 | p_g2) + plot_layout(widths = c(2.1, 0.8))))

  design <- "
AABB
CCDD
EEFF
GGGG
"
  main_figure <- panel_a + panel_b + panel_c + panel_d + panel_e + panel_f + panel_g +
    plot_layout(design = design, heights = c(0.9, 1.1, 1.1, 1.55)) +
    plot_annotation(
      title = "Transcriptomic reversal generates exploratory hypotheses\nwithout definitive cross-platform validation",
      subtitle = "Figure 8 v2 | EXTENDED_DATA_ONLY",
      caption = "All evidence is computational or cancer-cell based. No panel establishes treatment efficacy, normal-cell safety, or clinical actionability.",
      tag_levels = "a",
      theme = theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 7), plot.caption = element_text(size = 5.5))
    )
  figure8_v2_save_plot(main_figure, "figure8_v2_mainfigure_a_to_g", width_mm = 183, height_mm = 230, dpi = 600)
  saveRDS(main_figure, file.path(FIGURE8_V2_DATA, "figure8_v2_mainfigure_a_to_g_plot.rds"))
  invisible(main_figure)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  figure8_v2_plot_main()
  cat("FIGURE8_V2_MAIN_FIGURE exported\n")
}
