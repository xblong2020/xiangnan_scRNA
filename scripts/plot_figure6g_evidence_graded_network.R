#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages({library(igraph); library(tidygraph); library(ggraph)})

edges <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6g_edge_evidence.tsv"))
nodes <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6g_node_attributes.tsv"))[included_main_network == TRUE]
nodes[, plot_label := label]
nodes[node == "Stress-transition axis", plot_label := "Stress\ntransition\n(B)"]
nodes[node == "Proliferation", plot_label := "Prolif.\noutput"]
nodes[node == "CNV-associated\nmalignant signature", plot_label := "CNV\nsignature"]
edges <- edges[source %in% nodes$node & target %in% nodes$node]
set.seed(20260805)
graph <- tbl_graph(nodes = nodes, edges = edges[, .(from = source, to = target, evidence_score, evidence_grade, effect_direction)], directed = TRUE)
edge_lty <- c(strong = "solid", moderate = "solid", weak = "dashed", unresolved = "dotted", opposite = "dotdash")
edge_col <- c(positive = lancet_palette[2], negative = lancet_palette[1])
node_cols <- c(identity_axis = unname(axis_palette["identity_axis"]), stress_axis = unname(axis_palette["stress_axis"]),
  sox4_axis = unname(axis_palette["sox4_axis"]), outcome = unname(auxiliary_purple))
p <- ggraph(graph, layout = "circle") +
  geom_edge_link(aes(filter = evidence_grade != "unresolved", width = pmax(.25, evidence_score),
    linetype = evidence_grade, colour = effect_direction),
    arrow = grid::arrow(length = grid::unit(.075, "inches"), type = "closed"), end_cap = circle(7, "mm"), alpha = .9) +
  geom_edge_link(aes(filter = evidence_grade == "unresolved", width = pmax(.25, evidence_score), linetype = evidence_grade),
    colour = neutral_gray, alpha = .55) +
  geom_node_point(aes(fill = ifelse(is.na(axis), "outcome", axis)), shape = 21, size = 12.5, colour = dark_text, stroke = .6) +
  geom_node_text(aes(label = plot_label), size = 2, colour = "white", lineheight = .82) +
  scale_edge_width(range = c(.35, 2.1), name = "Evidence score") + scale_edge_linetype_manual(values = edge_lty, name = "Evidence grade") +
  scale_edge_colour_manual(values = edge_col, name = "Predicted sign") + scale_fill_manual(values = node_cols, guide = "none") +
  labs(title = "Evidence-graded directional network",
    caption = "Edge direction represents computational support and does not establish direct causality.") +
  coord_cartesian(xlim = c(-1.24, 1.24), ylim = c(-1.14, 1.14), clip = "off") +
  theme_void(base_family = "sans") + theme(plot.title = element_text(size=10,hjust=.5,colour=dark_text), plot.caption=element_text(size=7,hjust=.5,colour=dark_text),
    legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size=6.7), legend.title = element_text(size=7),
    legend.key.width = grid::unit(.65, "cm"), plot.margin=margin(8,10,8,10)) +
  guides(edge_linetype = guide_legend(nrow = 1), edge_width = guide_legend(nrow = 1), edge_colour = guide_legend(nrow = 1))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6g_evidence_graded_network")
figure6_save(p, out_dir, "figure6g_evidence_graded_network", 8.2, 5.6)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6g_plot.rds"))
