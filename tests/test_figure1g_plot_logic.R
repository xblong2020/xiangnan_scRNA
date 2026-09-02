helper_path <- file.path("scripts", "figure1g_stemness_helpers.R")
stopifnot(file.exists(helper_path))
source(helper_path)

state_order <- c("normal_hepatocyte_like", "stressed_injured_hepatocyte")
scores <- data.frame(
  cell_id = c("cell_1", "cell_2", "cell_3", "cell_4"),
  hepatocyte_state_label = c(
    "stressed_injured_hepatocyte",
    "normal_hepatocyte_like",
    "unknown_state",
    "normal_hepatocyte_like"
  ),
  CytoTRACE2_Score = c(0.20, 0.80, 0.50, NA_real_),
  stringsAsFactors = FALSE
)

prepared <- prepare_figure1g_data(scores, state_order)
stopifnot(nrow(prepared) == 2L)
stopifnot(identical(levels(prepared$hepatocyte_state_label), state_order))
stopifnot(identical(as.character(prepared$hepatocyte_state_label), c(
  "stressed_injured_hepatocyte", "normal_hepatocyte_like"
)))

summary_df <- summarize_figure1g_stemness(prepared, state_order)
stopifnot(identical(summary_df$n_cells, c(1L, 1L)))
stopifnot(isTRUE(all.equal(summary_df$median_cytotrace2, c(0.80, 0.20))))

log_scaled <- log10_minmax_scores(c(0, 0.09, 0.99), pseudocount = 0.01)
stopifnot(isTRUE(all.equal(range(log_scaled), c(0, 1))))
stopifnot(log_scaled[[2]] > log_scaled[[1]])
stopifnot(log_scaled[[3]] > log_scaled[[2]])

message("PASS: Figure 1G stemness data preparation logic")
