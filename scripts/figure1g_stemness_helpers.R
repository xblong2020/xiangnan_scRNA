prepare_figure1g_data <- function(scores, state_order) {
  required_columns <- c("cell_id", "hepatocyte_state_label", "CytoTRACE2_Score")
  missing_columns <- setdiff(required_columns, names(scores))
  if (length(missing_columns) > 0L) {
    stop("Scores table is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  out <- scores[, required_columns, drop = FALSE]
  out$CytoTRACE2_Score <- as.numeric(out$CytoTRACE2_Score)
  out <- out[
    is.finite(out$CytoTRACE2_Score) & out$hepatocyte_state_label %in% state_order,
    ,
    drop = FALSE
  ]
  out$hepatocyte_state_label <- factor(out$hepatocyte_state_label, levels = state_order)
  out
}

summarize_figure1g_stemness <- function(plot_data, state_order) {
  split_scores <- split(plot_data$CytoTRACE2_Score, plot_data$hepatocyte_state_label, drop = FALSE)
  data.frame(
    hepatocyte_state_label = factor(state_order, levels = state_order),
    n_cells = as.integer(vapply(split_scores, length, integer(1))),
    mean_cytotrace2 = as.numeric(vapply(split_scores, function(x) if (length(x) > 0L) mean(x) else NA_real_, numeric(1))),
    median_cytotrace2 = as.numeric(vapply(split_scores, function(x) if (length(x) > 0L) median(x) else NA_real_, numeric(1))),
    q1_cytotrace2 = as.numeric(vapply(split_scores, function(x) if (length(x) > 0L) unname(stats::quantile(x, 0.25)) else NA_real_, numeric(1))),
    q3_cytotrace2 = as.numeric(vapply(split_scores, function(x) if (length(x) > 0L) unname(stats::quantile(x, 0.75)) else NA_real_, numeric(1))),
    stringsAsFactors = FALSE
  )
}

log10_minmax_scores <- function(scores, pseudocount = 1e-4) {
  if (!is.numeric(scores) || !is.finite(pseudocount) || pseudocount <= 0) {
    stop("scores must be numeric and pseudocount must be a positive finite number.", call. = FALSE)
  }
  transformed <- log10(scores + pseudocount)
  if (any(!is.finite(transformed))) {
    stop("CytoTRACE2 scores must be non-negative finite values.", call. = FALSE)
  }
  value_range <- range(transformed)
  if (diff(value_range) == 0) return(rep(0, length(transformed)))
  (transformed - value_range[[1]]) / diff(value_range)
}
