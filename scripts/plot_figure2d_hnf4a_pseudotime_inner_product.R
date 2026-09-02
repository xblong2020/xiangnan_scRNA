#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))
root <- figure2_project_root()
c_dir <- figure2_get_arg("--figure2c-data-dir", file.path(root, "metadata/driver/figure2c_hnf4a"))
out_dir <- figure2_get_arg("--out-dir", file.path(root, "metadata/driver/figure2d_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2d_hnf4a"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")
k <- as.integer(figure2_get_arg("--k-neighbors", "50"))
n_bins <- as.integer(figure2_get_arg("--n-bins", "10"))
n_boot <- as.integer(figure2_get_arg("--n-bootstrap", "2000"))
seed <- as.integer(figure2_get_arg("--seed", "15071990"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim(gzfile(file.path(c_dir, "figure2c_hnf4a_cells_with_tsne_projection.tsv.gz")),
                    stringsAsFactors = FALSE)
ug <- read.delim(gzfile(file.path(c_dir, "figure2c_hnf4a_inner_product_grid_umap.tsv.gz")),
                 stringsAsFactors = FALSE)
tg <- read.delim(gzfile(file.path(c_dir, "figure2c_hnf4a_inner_product_grid_tsne.tsv.gz")),
                 stringsAsFactors = FALSE)

map_pt <- function(df, grid, x_col, y_col, k) {
  out <- grid; out$pseudotime_grid <- out$pseudotime_sd <- NA_real_
  for (i in seq_len(nrow(out))) {
    d <- sqrt((df[[x_col]] - out$grid_x[i])^2 + (df[[y_col]] - out$grid_y[i])^2)
    ix <- order(d)[seq_len(min(k, nrow(df)))]
    bw <- max(median(d[ix]), 1e-8)
    w <- exp(-(d[ix]^2) / (2 * bw^2)); w <- w / sum(w)
    mu <- sum(df$pseudotime[ix] * w)
    out$pseudotime_grid[i] <- mu
    out$pseudotime_sd[i] <- sqrt(sum(w * (df$pseudotime[ix] - mu)^2))
  }
  out$valid <- (out$keep_score %in% TRUE) & is.finite(out$inner_product_score_grid) &
    is.finite(out$pseudotime_grid)
  out
}
umap <- map_pt(cells, ug, "umap_1", "umap_2", k)
tsne <- map_pt(cells, tg, "tsne_1", "tsne_2", k)
write.table(umap, gzfile(file.path(out_dir, "figure2d_hnf4a_pseudotime_inner_product_umap.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tsne, gzfile(file.path(out_dir, "figure2d_hnf4a_pseudotime_inner_product_tsne.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summarize_bins <- function(df, space) {
  sub <- df[df$valid, , drop = FALSE]
  sub$pseudotime_bin <- cut(sub$pseudotime_grid, breaks = seq(0, 1, length.out = n_bins + 1),
                            include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(split(sub, sub$pseudotime_bin), function(x) data.frame(
    space = space, pseudotime_bin = unique(x$pseudotime_bin),
    pseudotime_center = mean(x$pseudotime_grid), n_grid_points = nrow(x),
    score_mean = mean(x$inner_product_score_grid),
    score_median = median(x$inner_product_score_grid),
    score_q25 = unname(quantile(x$inner_product_score_grid, .25)),
    score_q75 = unname(quantile(x$inner_product_score_grid, .75)),
    positive_fraction = mean(x$inner_product_score_grid > 0)
  )))
}
bins <- rbind(summarize_bins(umap, "CellOracle UMAP grid"),
              summarize_bins(tsne, "Expression t-SNE grid"))
write.table(bins, file.path(out_dir, "figure2d_hnf4a_pseudotime_bin_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bootstrap_rho <- function(df) {
  set.seed(seed)
  sub <- df[df$valid, , drop = FALSE]
  rho <- cor(sub$pseudotime_grid, sub$inner_product_score_grid, method = "spearman")
  boot <- replicate(n_boot, {
    ix <- sample.int(nrow(sub), nrow(sub), replace = TRUE)
    suppressWarnings(cor(sub$pseudotime_grid[ix], sub$inner_product_score_grid[ix],
                         method = "spearman"))
  })
  c(rho = unname(rho), ci_low = unname(quantile(boot, .025, na.rm = TRUE)),
    ci_high = unname(quantile(boot, .975, na.rm = TRUE)))
}

stage_stats <- function(df, space) {
  sub <- df[df$valid, , drop = FALSE]
  sub$stage <- cut(sub$pseudotime_grid, c(-Inf, 1/3, 2/3, Inf),
                   labels = c("early", "middle", "late"))
  stage <- do.call(rbind, lapply(split(sub, sub$stage), function(x) data.frame(
    space = space, stage = as.character(unique(x$stage)), n_grid_points = nrow(x),
    score_mean = mean(x$inner_product_score_grid), score_median = median(x$inner_product_score_grid),
    score_q25 = unname(quantile(x$inner_product_score_grid, .25)),
    score_q75 = unname(quantile(x$inner_product_score_grid, .75)),
    positive_fraction = mean(x$inner_product_score_grid > 0)
  )))
  kw <- kruskal.test(inner_product_score_grid ~ stage, data = sub)
  jt <- suppressWarnings(cor.test(as.numeric(sub$stage), sub$inner_product_score_grid,
                                  method = "spearman", exact = FALSE))
  list(table = stage, kruskal_p = unname(kw$p.value),
       ordered_trend_rho = unname(jt$estimate), ordered_trend_p = unname(jt$p.value))
}
u_stage <- stage_stats(umap, "CellOracle UMAP grid")
t_stage <- stage_stats(tsne, "Expression t-SNE grid")
write.table(rbind(u_stage$table, t_stage$table),
            file.path(out_dir, "figure2d_hnf4a_stage_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_relation <- function(df, bin_df, title, tag = NULL) {
  sub <- df[df$valid, , drop = FALSE]
  lim <- max(abs(sub$inner_product_score_grid), na.rm = TRUE)
  ggplot(sub, aes(pseudotime_grid, inner_product_score_grid)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .35, colour = "grey45") +
    geom_point(aes(colour = inner_product_score_grid), size = 1.45, alpha = .78) +
    geom_line(data = bin_df[order(bin_df$pseudotime_center), ],
              aes(pseudotime_center, score_median), inherit.aes = FALSE,
              colour = "black", linewidth = .55) +
    geom_point(data = bin_df, aes(pseudotime_center, score_median), inherit.aes = FALSE,
               shape = 21, size = 2, stroke = .45, fill = "white", colour = "black") +
    scale_colour_gradient2(low = figure2_diverging["low"], mid = figure2_diverging["mid"],
      high = figure2_diverging["high"], midpoint = 0, limits = c(-lim, lim),
      oob = scales::squish, name = "Inner product score (PS)") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    labs(x = "driver_main_strict__pseudotime_rank", y = "Inner product score",
         title = title, tag = tag) + figure2_theme()
}
ub <- bins[bins$space == "CellOracle UMAP grid", ]
tb <- bins[bins$space == "Expression t-SNE grid", ]
pu <- plot_relation(umap, ub, "HNF4A perturbation score along pseudotime", "Figure 2D")
pt <- plot_relation(tsne, tb, "HNF4A perturbation score along pseudotime (t-SNE sensitivity)")
figure2_save(pu, figure_dir, "figure2d_hnf4a_pseudotime_inner_product_umap", 5.7, 4.3)
figure2_save(pt, figure_dir, "figure2d_hnf4a_pseudotime_inner_product_tsne", 5.7, 4.3)
figure2_save(pu + pt + plot_layout(guides = "collect") & theme(legend.position = "right"),
             figure_dir, "figure2d_hnf4a_pseudotime_inner_product_umap_tsne", 10.8, 4.3)

urho <- bootstrap_rho(umap); trho <- bootstrap_rho(tsne)
report <- list(
  module = "Figure 2D", target_tf = target_tf,
  analysis = "HNF4A perturbation score along pseudotime",
  score_definition = "HNF4A perturbation vector dot baseline developmental vector",
  n_cells = nrow(cells), n_bins = n_bins, bootstrap_replicates = n_boot, seed = seed,
  umap = list(n_grid_points = sum(umap$valid), score_range = range(umap$inner_product_score_grid[umap$valid]),
              spearman_rho = urho["rho"], spearman_ci95 = urho[c("ci_low", "ci_high")],
              early_middle_late_kruskal_p = u_stage$kruskal_p,
              early_to_late_ordered_trend_rho = u_stage$ordered_trend_rho,
              early_to_late_ordered_trend_p = u_stage$ordered_trend_p),
  tsne = list(n_grid_points = sum(tsne$valid), score_range = range(tsne$inner_product_score_grid[tsne$valid]),
              spearman_rho = trho["rho"], spearman_ci95 = trho[c("ci_low", "ci_high")],
              early_middle_late_kruskal_p = t_stage$kruskal_p,
              early_to_late_ordered_trend_rho = t_stage$ordered_trend_rho,
              early_to_late_ordered_trend_p = t_stage$ordered_trend_p),
  caveat = "Direction and timing are reported from the observed predicted perturbation without assuming a positive or early effect."
)
figure2_write_json(report, file.path(out_dir, "figure2d_hnf4a_report.json"))
