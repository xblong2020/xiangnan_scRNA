options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 1800)

message("R version: ", getRversion())
message("CRAN repo: ", getOption("repos")["CRAN"])

cran_pkgs <- c("BH", "TH.data", "multcomp")
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
}

message("Installing infercnv")
BiocManager::install("infercnv", ask = FALSE, update = FALSE)

pkgs <- c("BH", "TH.data", "multcomp", "infercnv")
message("Final retry package status:")
for (pkg in pkgs) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(packageVersion(pkg)) else ""
  message(pkg, "\t", available, "\t", version)
}

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("infercnv installation did not complete.")
}
