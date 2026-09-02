#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages(library(patchwork))

fit <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6f_model_fit_summary.tsv"))
model_cols <- c("Model 1: Linear cascade" = lancet_palette[1], "Model 2: Parallel convergence" = lancet_palette[5], "Model 3: Partially ordered" = lancet_palette[3])
node_labels <- c(A="Identity\nloss", B="Stress\ntransition", C="SOX4\nstabilization", F="Malignant\nfate")
nd <- rbindlist(list(
  data.table(model=names(model_cols)[1], id=c("A","B","C","F"), x=1:4, y=1),
  data.table(model=names(model_cols)[2], id=c("A","B","C","F"), x=c(1,1,1,3.6), y=c(1.38,1,.62,1)),
  data.table(model=names(model_cols)[3], id=c("A","B","C","F"), x=1:4, y=c(1.25,.78,1.25,.78))
))
nd[, node := node_labels[id]]
edges <- rbindlist(list(
  data.table(model=names(model_cols)[1], from=c("A","B","C"), to=c("B","C","F")),
  data.table(model=names(model_cols)[2], from=c("A","B","C"), to=c("F","F","F")),
  data.table(model=names(model_cols)[3], from=c("A","A","B","A","B","C"), to=c("B","C","C","F","F","F"))
))
edges <- merge(edges, nd[, .(model, from=id, x_from=x, y_from=y)], by=c("model","from"))
edges <- merge(edges, nd[, .(model, to=id, x_to=x, y_to=y)], by=c("model","to"))
edges[, edge_distance := sqrt((x_to-x_from)^2 + (y_to-y_from)^2)]
edges[, `:=`(
  x_start = x_from + .24 * (x_to-x_from)/edge_distance,
  y_start = y_from + .24 * (y_to-y_from)/edge_distance,
  x_end = x_to - .30 * (x_to-x_from)/edge_distance,
  y_end = y_to - .30 * (y_to-y_from)/edge_distance
)]
p_net <- ggplot() +
  geom_segment(data = edges, aes(x_start, y_start, xend = x_end, yend = y_end), colour = neutral_gray,
    arrow = grid::arrow(length = grid::unit(.045, "inches"), type = "closed"), linewidth = .45) +
  geom_label(data = nd, aes(x, y, label = node, fill = model), size = 2.3, label.size = .25, colour = dark_text) +
  facet_wrap(~model, ncol = 1) + scale_fill_manual(values = model_cols, guide = "none") +
  coord_cartesian(xlim = c(.55,4.45), ylim = c(.42,1.58), clip = "off") +
  theme_void(base_family = "sans") + theme(strip.text = element_text(size = 7.5, colour = dark_text), plot.title = element_text(size=9,hjust=.5)) +
  labs(title = "Predefined identifiable models")
metrics <- melt(fit, id.vars = c("model", "selected"), measure.vars = c("aic", "bic", "cv_rmse", "mean_edge_stability"),
  variable.name = "metric", value.name = "value")
metrics[, performance := ifelse(metric == "mean_edge_stability", figure6_safe_scale(value), -figure6_safe_scale(value))]
metrics[, metric := factor(metric, levels = c("aic","bic","cv_rmse","mean_edge_stability"), labels = c("AIC","BIC","CV RMSE","Edge stability"))]
p_perf <- ggplot(metrics, aes(performance, metric, colour = model, shape = selected, group = model)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .35, colour = neutral_gray) +
  geom_point(size = 2.8, position = position_dodge(width = .48)) +
  scale_colour_manual(values = model_cols, name = NULL) + scale_shape_manual(values = c(`FALSE`=16, `TRUE`=8), name = "Selected by\nrank-sum rule") +
  labs(title = "Multi-metric model comparison", x = "Normalized performance (higher is better)", y = NULL,
    caption = "Selection: minimum rank sum of AIC, BIC, repeated-CV RMSE and bootstrap edge stability.") +
  figure6_theme() + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2))
p <- p_net + p_perf + plot_layout(widths = c(1.1, 1.45)) + plot_annotation(title = "Competing network models", theme = theme(plot.title=element_text(size=10,hjust=.5)))
out_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6f_competing_network_models")
figure6_save(p, out_dir, "figure6f_competing_network_models", 10.2, 5.2)
saveRDS(p, file.path(FIGURE6_METADATA_DIR, "figure6f_plot.rds"))
