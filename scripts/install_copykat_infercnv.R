options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = max(600, getOption("timeout", 60)))

message("R version: ", getRversion())
message("Library paths:")
message(paste(.libPaths(), collapse = "\n"))

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
}

install_cran_if_missing(c("BiocManager", "remotes"))

message("BiocManager version: ", packageVersion("BiocManager"))
message("Bioconductor version: ", BiocManager::version())

bioc_pkgs <- c("DNAcopy", "infercnv")
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

install_copykat <- function() {
  if (requireNamespace("copykat", quietly = TRUE)) {
    return(invisible(TRUE))
  }
  message("Installing copykat from GitHub via remotes::install_github")
  ok <- tryCatch({
    remotes::install_github("navinlabcode/copykat", dependencies = TRUE, upgrade = "never")
    TRUE
  }, error = function(e) {
    message("remotes::install_github failed: ", conditionMessage(e))
    FALSE
  })
  if (ok && requireNamespace("copykat", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Installing copykat from GitHub codeload zip fallback")
  tmp <- file.path(tempdir(), "copykat-main.zip")
  urls <- c(
    "https://github.com/navinlabcode/copykat/archive/refs/heads/main.zip",
    "https://github.com/navinlabcode/copykat/archive/refs/heads/master.zip"
  )
  downloaded <- FALSE
  for (url in urls) {
    message("Trying: ", url)
    downloaded <- tryCatch({
      download.file(url, tmp, mode = "wb", quiet = FALSE)
      file.exists(tmp) && file.info(tmp)$size > 0
    }, error = function(e) {
      message("download failed: ", conditionMessage(e))
      FALSE
    })
    if (downloaded) break
  }
  if (!downloaded) {
    stop("Could not download copykat zip from GitHub codeload.")
  }
  remotes::install_local(tmp, dependencies = TRUE, upgrade = "never")
  invisible(requireNamespace("copykat", quietly = TRUE))
}

install_copykat()

pkgs <- c("BiocManager", "remotes", "DNAcopy", "infercnv", "copykat")
message("Final package status:")
for (pkg in pkgs) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(packageVersion(pkg)) else ""
  message(pkg, "\t", available, "\t", version)
}

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("infercnv installation did not complete.")
}
if (!requireNamespace("copykat", quietly = TRUE)) {
  stop("copykat installation did not complete.")
}
