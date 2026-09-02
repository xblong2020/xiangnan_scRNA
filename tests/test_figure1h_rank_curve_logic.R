helper_path <- file.path("scripts", "figure1h_rank_curve_helpers.R")
stopifnot(file.exists(helper_path))
source(helper_path)

rss <- data.frame(
  state = rep(c("state_a", "state_b"), each = 6),
  regulon = rep(paste0("TF", 1:6, "(+)"), 2),
  rss = c(0.2, 0.6, 0.4, 0.5, 0.3, 0.1, 0.7, 0.2, 0.6, 0.1, 0.5, 0.4),
  stringsAsFactors = FALSE
)

ranked <- prepare_rss_rank_curves(rss, state_order = c("state_a", "state_b"), top_n = 5)
stopifnot(nrow(ranked) == 12L)
stopifnot(all(tapply(ranked$rank, ranked$state, min) == 1L))
stopifnot(all(tapply(!is.na(ranked$top_label), ranked$state, sum) == 5L))
stopifnot(ranked$regulon[ranked$state == "state_a" & ranked$rank == 1L] == "TF2(+)")

labels <- build_top_label_positions(ranked)
stopifnot(nrow(labels) == 10L)
stopifnot(all(tapply(labels$label_x, labels$state, function(x) length(unique(x))) == 5L))
stopifnot(all(tapply(labels$label_y, labels$state, function(x) length(unique(x))) == 5L))
stopifnot(all(labels$point_x == labels$rank))
stopifnot(all(labels$point_y == labels$rss))
stopifnot(all(tapply(labels$label_y, labels$state, function(x) all(diff(x) < 0))))
stopifnot(identical(lancet_state_palette(6), c(
  "#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF"
)))

message("PASS: Figure 1H RSS rank-curve logic")
