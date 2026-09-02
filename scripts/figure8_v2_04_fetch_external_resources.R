#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

figure8_v2_curl_handle <- function() {
  if (!requireNamespace("curl", quietly = TRUE)) stop("R package 'curl' is required for official external resources")
  handle <- curl::new_handle()
  curl::handle_setopt(handle, useragent = "Figure8V2-HCC-reproducible-research/1.0")
  handle
}

figure8_v2_powershell_fetch <- function(url, destination, json = FALSE) {
  powershell <- Sys.which("powershell.exe")
  if (!nzchar(powershell)) powershell <- Sys.which("powershell")
  if (!nzchar(powershell)) stop("PowerShell is required for Figshare downloads in this Windows environment")
  script <- tempfile(fileext = ".ps1")
  on.exit(unlink(script), add = TRUE)
  lines <- if (json) c(
    "param([string]$Url,[string]$OutFile)",
    "$ErrorActionPreference='Stop'", "$ProgressPreference='SilentlyContinue'",
    "Invoke-RestMethod -Uri $Url | ConvertTo-Json -Depth 50 | Set-Content -LiteralPath $OutFile -Encoding utf8"
  ) else c(
    "param([string]$Url,[string]$OutFile)",
    "$ErrorActionPreference='Stop'", "$ProgressPreference='SilentlyContinue'",
    "Invoke-WebRequest -Uri $Url -OutFile $OutFile"
  )
  writeLines(lines, script, useBytes = TRUE)
  output <- system2(
    powershell,
    args = c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", shQuote(script), "-Url", shQuote(url), "-OutFile", shQuote(destination)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status") %||% 0L
  if (status != 0L || !file.exists(destination)) stop("PowerShell fetch failed: ", paste(output, collapse = "\n"))
  invisible(destination)
}

figure8_v2_fetch_json <- function(url) {
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  figure8_v2_powershell_fetch(url, path, json = TRUE)
  jsonlite::fromJSON(path)
}

figure8_v2_external_resource_plan <- function() {
  rbindlist(list(
    data.table(
      resource_id = "prism_23q2",
      article_id = 23600310L,
      release = "Repurposing Public 23Q2 v4",
      doi = "10.6084/m9.figshare.23600310.v4",
      api_url = "https://api.figshare.com/v2/articles/23600310",
      filename = c(
        "Repurposing_Public_23Q2_Cell_Line_Meta_Data.csv",
        "Repurposing_Public_23Q2_Readme.txt",
        "Repurposing_Public_23Q2_Treatment_Meta_Data.csv",
        "Repurposing_Public_23Q2_LFC_COLLAPSED.csv",
        "Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv",
        "Repurposing_Public_23Q2_QC_table.csv"
      ),
      role = c("cell_line_metadata", "readme", "treatment_metadata", "primary_single_dose_lfc", "compound_annotation", "quality_control")
    ),
    data.table(
      resource_id = "prism_19q4_secondary",
      article_id = 9393293L,
      release = "PRISM Repurposing 19Q4 v4 secondary",
      doi = "10.6084/m9.figshare.9393293.v4",
      api_url = "https://api.figshare.com/v2/articles/9393293",
      filename = c(
        "secondary-screen-readme.txt",
        "secondary-screen-dose-response-curve-parameters.csv",
        "secondary-screen-replicate-collapsed-logfold-change.csv",
        "secondary-screen-replicate-collapsed-treatment-info.csv",
        "secondary-screen-cell-line-info.csv"
      ),
      role = c("readme", "secondary_dose_response", "secondary_collapsed_lfc", "secondary_treatment_metadata", "secondary_cell_line_metadata")
    ),
    data.table(
      resource_id = "depmap_23q2_model_metadata",
      article_id = 22765112L,
      release = "DepMap 23Q2 Public v4",
      doi = "10.6084/m9.figshare.22765112.v4",
      api_url = "https://api.figshare.com/v2/articles/22765112",
      filename = "Model.csv",
      role = "release_matched_model_tissue_disease_metadata"
    )
  ), use.names = TRUE)
}

figure8_v2_plan_for_article <- function(plan, requested_article_id) {
  as.data.table(plan)[article_id == as.integer(requested_article_id)]
}

figure8_v2_validate_external_file <- function(path, expected_size, expected_md5) {
  if (!file.exists(path)) return(data.table(status = "missing", observed_size = NA_real_, observed_md5 = NA_character_))
  observed_size <- as.numeric(file.info(path)$size)
  if (is.finite(expected_size) && observed_size != as.numeric(expected_size)) {
    return(data.table(status = "size_mismatch", observed_size = observed_size, observed_md5 = NA_character_))
  }
  observed_md5 <- unname(tools::md5sum(path))
  status <- if (!is.na(expected_md5) && nzchar(expected_md5) && observed_md5 != expected_md5) "md5_mismatch" else "verified"
  data.table(status = status, observed_size = observed_size, observed_md5 = observed_md5)
}

figure8_v2_missing_download_rows <- function(result) {
  result <- as.data.table(copy(result))
  if (!"error" %in% names(result)) result[, error := NA_character_]
  result[status != "verified", .(resource_id, release, filename, role, status, error)]
}

figure8_v2_download_one <- function(url, destination, expected_size, expected_md5, retries = 3L) {
  if (file.exists(destination)) {
    existing <- figure8_v2_validate_external_file(destination, expected_size, expected_md5)
    if (existing$status == "verified") return(existing)
    stop("Existing external file fails checksum and will not be overwritten automatically: ", destination)
  }
  partial <- paste0(destination, ".partial")
  on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    result <- tryCatch({
      figure8_v2_powershell_fetch(url, partial, json = FALSE)
      check <- figure8_v2_validate_external_file(partial, expected_size, expected_md5)
      if (check$status != "verified") stop("Downloaded file verification failed: ", check$status)
      if (!file.rename(partial, destination)) stop("Could not atomically rename completed download")
      figure8_v2_validate_external_file(destination, expected_size, expected_md5)
    }, error = function(e) e)
    if (!inherits(result, "error")) return(result)
    last_error <- conditionMessage(result)
    if (file.exists(partial)) unlink(partial)
  }
  data.table(status = "download_failed", observed_size = NA_real_, observed_md5 = NA_character_, error = last_error)
}

figure8_v2_fetch_external_resources <- function() {
  figure8_v2_init_dirs()
  raw_dir <- file.path(FIGURE8_V2_DATA, "external_raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  plan <- figure8_v2_external_resource_plan()
  article_rows <- list()
  for (article_id in unique(plan$article_id)) {
    article_plan <- figure8_v2_plan_for_article(plan, article_id)
    api_url <- unique(article_plan$api_url)
    article <- figure8_v2_fetch_json(api_url)
    files <- as.data.table(article$files)
    selected <- merge(article_plan, files[, .(
      filename = name, figshare_file_id = id, expected_size = size,
      download_url, expected_md5 = supplied_md5, mimetype
    )], by = "filename", all.x = TRUE)
    selected[, `:=`(
      article_title = article$title,
      article_version = article$version,
      article_published_date = article$published_date,
      article_modified_date = article$modified_date,
      license = article$license$name
    )]
    article_rows[[as.character(article_id)]] <- selected
    write_json(article, file.path(raw_dir, paste0("figshare_article_", article_id, ".json")), pretty = TRUE, auto_unbox = TRUE)
  }
  manifest <- rbindlist(article_rows, fill = TRUE)
  manifest[, local_path := file.path(raw_dir, filename)]
  result_rows <- vector("list", nrow(manifest))
  for (idx in seq_len(nrow(manifest))) {
    row <- manifest[idx]
    if (is.na(row$download_url) || !nzchar(row$download_url)) {
      check <- data.table(status = "metadata_missing_download_url", observed_size = NA_real_, observed_md5 = NA_character_, error = "Figshare file not found in article metadata")
    } else {
      check <- figure8_v2_download_one(row$download_url, row$local_path, row$expected_size, row$expected_md5)
    }
    result_rows[[idx]] <- cbind(row, check)
  }
  result <- rbindlist(result_rows, fill = TRUE)
  result[, retrieval_date_utc := format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]

  cellosaurus_url <- "https://api.cellosaurus.org/release-info?format=json"
  cellosaurus_path <- file.path(raw_dir, "cellosaurus_release_info.json")
  cellosaurus_status <- tryCatch({
    release <- figure8_v2_fetch_json(cellosaurus_url)
    write_json(release, cellosaurus_path, pretty = TRUE, auto_unbox = TRUE)
    "verified"
  }, error = function(e) paste0("download_failed: ", conditionMessage(e)))

  figure8_v2_write_tsv(result, "figure8_v2_external_resource_manifest.tsv")
  missing <- figure8_v2_missing_download_rows(result)
  missing <- rbindlist(list(
    missing,
    data.table(
      resource_id = "repurposing_hub_standalone",
      release = "CLUE Repurposing Hub 2018",
      filename = NA_character_, role = "curated_moa_target_annotation",
      status = "standalone_download_not_available",
      error = "CLUE site retired in 2026; PRISM treatment/compound metadata and official exact APIs are used instead"
    ),
    data.table(
      resource_id = "cellosaurus",
      release = "current API release info",
      filename = basename(cellosaurus_path), role = "cell_line_verification",
      status = if (cellosaurus_status == "verified") "verified" else "download_failed",
      error = if (cellosaurus_status == "verified") NA_character_ else cellosaurus_status
    )
  ), fill = TRUE)
  figure8_v2_write_tsv(missing, "figure8_v2_missing_external_resource_audit.tsv")
  if (any(result$status %in% c("size_mismatch", "md5_mismatch"))) stop("External resource checksum validation failed")
  invisible(result)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_fetch_external_resources()
  cat("FIGURE8_V2_EXTERNAL verified=", sum(result$status == "verified"), "/", nrow(result), "\n", sep = "")
}
