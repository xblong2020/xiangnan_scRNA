prepare_rss_rank_curves <- function(rss, state_order, top_n = 5L) {
  required <- c("state", "regulon", "rss")
  missing <- setdiff(required, names(rss))
  if (length(missing) > 0L) {
    stop("RSS table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  absent_states <- setdiff(state_order, unique(as.character(rss$state)))
  if (length(absent_states) > 0L) {
    stop("RSS table is missing states: ", paste(absent_states, collapse = ", "), call. = FALSE)
  }

  pieces <- lapply(state_order, function(state_name) {
    state_df <- rss[as.character(rss$state) == state_name, , drop = FALSE]
    state_df$rss <- as.numeric(state_df$rss)
    state_df <- state_df[order(state_df$rss, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    state_df$rank <- seq_len(nrow(state_df))
    state_df$top_label <- ifelse(state_df$rank <= as.integer(top_n), state_df$regulon, NA_character_)
    state_df
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out$state <- factor(out$state, levels = state_order)
  out
}

build_top_label_positions <- function(ranked, top_pad = 0.18, bottom_pad = 0.02) {
  top <- ranked[!is.na(ranked$top_label), , drop = FALSE]
  if (nrow(top) == 0L) {
    stop("No labelled regulons were found.", call. = FALSE)
  }

  state_levels <- if (is.factor(ranked$state)) levels(ranked$state) else unique(as.character(ranked$state))
  pieces <- lapply(state_levels, function(state_name) {
    state_all <- ranked[as.character(ranked$state) == state_name, , drop = FALSE]
    state_top <- top[as.character(top$state) == state_name, , drop = FALSE]
    if (nrow(state_top) == 0L) return(NULL)

    state_top <- state_top[order(state_top$rank), , drop = FALSE]
    rss_span <- diff(range(state_all$rss, na.rm = TRUE))
    if (!is.finite(rss_span) || rss_span == 0) {
      rss_span <- max(abs(state_all$rss), na.rm = TRUE) * 0.1
    }
    if (!is.finite(rss_span) || rss_span == 0) rss_span <- 0.1

    state_top$point_x <- state_top$rank
    state_top$point_y <- state_top$rss
    x_end <- max(1.5, min(max(state_all$rank, na.rm = TRUE) * 0.75, 26.0))
    state_top$label_x <- seq(1.5, x_end, length.out = nrow(state_top))
    state_top$label_y <- max(state_all$rss, na.rm = TRUE) +
      seq(top_pad, bottom_pad, length.out = nrow(state_top)) * rss_span
    state_top
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out$state <- factor(out$state, levels = state_levels)
  out
}

lancet_state_palette <- function(n) {
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop("Package 'ggsci' is required.", call. = FALSE)
  }
  ggsci::pal_lancet()(as.integer(n))
}
