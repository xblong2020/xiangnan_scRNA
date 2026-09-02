root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_01_build_continuous_signature.R")
if (!file.exists(script)) stop("Expected RED failure: continuous-signature module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

weights <- figure8_v2_component_weights()
stopifnot(isTRUE(all.equal(sum(weights), 1)))
stopifnot(identical(names(weights), c(
  "state_component", "trajectory_component", "axis_A_component", "axis_B_component",
  "axis_C_component", "malignant_state_component", "perturbation_component"
)))

synthetic <- data.table::data.table(
  gene = c("IDENTITY", "STRESS", "CONFLICT", "UNSUPPORTED"),
  state_component = c(0.8, -0.8, 0.8, NA),
  trajectory_component = c(0.6, -0.6, -0.8, NA),
  axis_A_component = c(1, 0, 1, 0),
  axis_B_component = c(0, -1, -1, 0),
  axis_C_component = c(0, -0.5, 0, 0),
  malignant_state_component = c(0, -0.8, 0, 0),
  perturbation_component = c(0.5, -0.5, 0, 0)
)
combined <- figure8_v2_combine_components(synthetic)
stopifnot(combined[gene == "IDENTITY", final_rescue_vscore] > 0)
stopifnot(combined[gene == "STRESS", final_rescue_vscore] < 0)
stopifnot(combined[gene == "CONFLICT", conflict_flag])
stopifnot(combined[gene == "UNSUPPORTED", final_rescue_vscore] == 0)
stopifnot(all(combined$final_rescue_vscore >= -1 & combined$final_rescue_vscore <= 1))

balance <- figure8_v2_axis_balance(data.table::data.table(
  gene = c("A1", "A2", "B1", "C1"),
  axis_A_component = c(1, 1, 0, 0),
  axis_B_component = c(0, 0, -1, 0),
  axis_C_component = c(0, 0, 0, -1)
))
stopifnot(isTRUE(all.equal(sum(balance$absolute_mass_fraction), 1)))
stopifnot(!any(balance$severe_axis_domination))

model <- figure8_v2_model_genes(root)
stopifnot(nrow(model) == 978L, model$gene[1] == "GNPDA1")
stopifnot(identical(model$model_gene_order, seq_len(978L)))

toy <- Matrix::Matrix(matrix(c(1, 0, 3, 2, 4, 0), nrow = 3, dimnames = list(NULL, c("GENE1", "GENE2"))), sparse = TRUE)
toy_means <- figure8_v2_sample_means(toy, c("S1", "S1", "S2"), minimum_cells = 1L)
stopifnot(all(c("sample_id", "n_cells", "GENE1", "GENE2") %in% names(toy_means)))
stopifnot(toy_means[sample_id == "S1", n_cells] == 2L)

toy_long <- figure8_v2_continuous_long(combined)
stopifnot("desired_direction" %in% names(toy_long))
stopifnot(!any(toy_long$desired_direction == "zero"))
stopifnot(!"final_rescue_direction" %in% names(toy_long))

cat("figure8_v2 signature logic tests passed\n")
