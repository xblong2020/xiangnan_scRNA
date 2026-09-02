#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

weights <- c(celloracle = .25, sctenifold = .15, temporal = .15, asymmetry = .15, model = .15, external = .10, control = .05)
effects <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_effects.tsv.gz"))[availability == "Available"]
asym <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry.tsv"))
sct <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "module9_2_network_directionality_matrix.tsv"))
order <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "module9_1_bootstrap_order_tests.tsv"))[run_id == "main_strict"]
fit <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6f_model_fit_summary.tsv"))
stab <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6f_bootstrap_edge_stability.tsv"))
selected <- fit[selected == TRUE, model][1]
conf_path <- file.path(FIGURE6_METADATA_DIR, "figure6_confounder_adjustment_summary.tsv")
conf <- if (file.exists(conf_path)) figure6_fread(conf_path) else data.table()

edges <- data.table(
  source = c("Identity axis","Identity axis","Stress-transition axis","Identity axis","Stress-transition axis","SOX4 axis",
    rep(c("Identity axis","Stress-transition axis","SOX4 axis"), 2)),
  target = c("Stress-transition axis","SOX4 axis","SOX4 axis","Malignant fate","Malignant fate","Malignant fate",
    rep("Proliferation",3), rep("CNV-associated\nmalignant signature",3)),
  source_code = c("A","A","B","A","B","C", "A","B","C", "A","B","C"),
  target_code = c("B","C","C","Fate","Fate","Fate", rep("Prolif",3), rep("CNV",3)),
  target_output = c("stress_transition_change","sox4_programme_change","sox4_programme_change",
    rep("malignant_fate_change",3), rep("proliferation_change",3), rep("cnv_malignant_signature_change",3))
)
axis_tfs <- list(A = c("HNF4A","PPARA"), B = c("EGR1","CEBPB","AP1_AGGREGATE"), C = "SOX4")
axis_tfs_sct <- list(A = c("HNF4A","PPARA"), B = unique(c("EGR1","CEBPB",FIGURE6_AP1_MEMBERS)), C = "SOX4")
axis_map_sct <- c(A = "A_upstream", B = "B_transition", C = "C_sox4")
comparison_map <- c("A_B"="Axis A → Axis B", "A_C"="Axis A → Axis C", "B_C"="Axis B → Axis C")
temporal_map <- c("A_B"="A_loss_before_B_transition", "B_C"="B_transition_before_C_sox4", "C_Fate"="C_sox4_before_or_equal_malignant_fate")

for (i in seq_len(nrow(edges))) {
  e <- edges[i]; key <- paste(e$source_code, e$target_code, sep = "_")
  eff <- effects[tf %in% axis_tfs[[e$source_code]] & output == e$target_output]
  controls <- effects[axis == "control" & output == e$target_output]
  candidate_mag <- mean(abs(eff$effect_estimate), na.rm = TRUE); control_mag <- median(abs(controls$effect_estimate), na.rm = TRUE)
  co_support <- if (!nrow(eff)) 0 else if (any(eff$fdr < .05, na.rm = TRUE) && candidate_mag > control_mag) 1 else if (candidate_mag > control_mag) .5 else .25
  sct_support <- 0
  if (e$target_code %in% c("A","B","C")) {
    q <- sct[perturb_tf %in% axis_tfs_sct[[e$source_code]] & target_axis == axis_map_sct[e$target_code]]
    sct_support <- if (nrow(q)) min(1, mean(q$sig_fraction, na.rm = TRUE) * 2) else 0
  }
  temporal_support <- 0
  if (key %in% names(temporal_map)) {
    q <- order[comparison == temporal_map[key]]
    temporal_support <- if (nrow(q)) mean(q$order_probability, na.rm = TRUE) else 0
  }
  asym_support <- 0
  if (key %in% names(comparison_map)) {
    q <- asym[comparison == comparison_map[key]]
    asym_support <- if (!nrow(q)) 0 else if (q$classification == "Forward-dominant") 1 else if (q$classification == "Symmetric/unresolved") .25 else 0
  }
  model_edge <- paste(e$source_code, "→", ifelse(e$target_code == "Fate", "fate", e$target_code))
  qst <- stab[model == selected & edge == model_edge]
  model_support <- if (nrow(qst)) qst$bootstrap_sign_stability[1] else 0
  control_support <- if (is.finite(candidate_mag) && is.finite(control_mag) && candidate_mag > control_mag) min(1, candidate_mag/(candidate_mag+control_mag+1e-9)) else 0
  qconf <- if (nrow(conf)) conf[tf %in% axis_tfs[[e$source_code]] & output == e$target_output & status == "estimable"] else data.table()
  adjustment_support <- if (nrow(qconf) && is.finite(mean(eff$effect_estimate, na.rm=TRUE)))
    mean(sign(qconf$adjusted_effect) == sign(mean(eff$effect_estimate, na.rm=TRUE)), na.rm=TRUE) else 0
  values <- c(celloracle = co_support, sctenifold = sct_support, temporal = temporal_support,
    asymmetry = asym_support, model = model_support, external = 0, control = control_support)
  score <- sum(weights * values)
  n_methods <- sum(values[c("celloracle","sctenifold","temporal","asymmetry","model")] >= .5)
  grade <- if (score >= .75 && n_methods >= 2 && temporal_support >= .5 && model_support >= .7 && asym_support >= .5 && adjustment_support >= .5) "strong"
    else if (score >= .50 && n_methods >= 2) "moderate" else if (score >= .30) "weak" else "unresolved"
  edges[i, names(values) := as.list(values)]
  edges[i, `:=`(evidence_score = score, n_supporting_methods = n_methods, evidence_grade = grade,
    adjustment_support = adjustment_support,
    effect_direction = ifelse(mean(eff$effect_estimate, na.rm = TRUE) >= 0, "positive", "negative"), selected_model = selected)]
}
figure6_fwrite(edges, file.path(FIGURE6_METADATA_DIR, "figure6g_edge_evidence.tsv"), compress = FALSE)
nodes <- data.table(
  node = c("Identity axis","Stress-transition axis","SOX4 axis","Malignant fate","Proliferation","CNV-associated\nmalignant signature"),
  node_type = c("axis","axis","axis","outcome","outcome","outcome"),
  axis = c("identity_axis","stress_axis","sox4_axis",NA,NA,NA),
  label = c("Identity\n(A)", "Stress transition\n(B)", "SOX4\n(C)", "Malignant\nfate", "Proliferation", "CNV-associated\nsignature"),
  included_main_network = TRUE
)
figure6_fwrite(nodes, file.path(FIGURE6_METADATA_DIR, "figure6g_node_attributes.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6G", weights = as.list(weights), selected_model = selected, layout_seed = 20260805,
  grade_rules = list(strong = "score>=0.75, >=2 methods, temporal>=0.5, model>=0.7, asymmetry>=0.5 and adjustment support>=0.5",
    moderate = "score>=0.50 and >=2 methods", weak = "score>=0.30", unresolved = "score<0.30 or insufficient support"),
  unavailable_edge_components = c("external edge-direction evidence", "perturb-malignant-fate reverse experiment", "restore/OE"),
  caveat = "Edge direction represents computational support and does not establish direct causality."
), file.path(FIGURE6_METADATA_DIR, "figure6g_network_report.json"))
