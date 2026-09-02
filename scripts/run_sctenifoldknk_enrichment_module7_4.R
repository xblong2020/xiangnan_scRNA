suppressPackageStartupMessages({
  library(dplyr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    perturbation = "metadata/driver/sctenifoldknk_module7_2_driver_union_all_perturbation_genes.tsv",
    background_genes = "data/processed/driver/sctenifoldknk_module7_1/driver_union_all/sctenifoldknk_genes.tsv",
    metadata_dir = "metadata/driver",
    fdr = 0.05,
    min_genes = 5
  )
  for (item in args) {
    parts <- strsplit(sub("^--", "", item), "=", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] %in% names(out)) {
      out[[parts[1]]] <- parts[2]
    }
  }
  out$fdr <- as.numeric(out$fdr)
  out$min_genes <- as.integer(out$min_genes)
  out
}

first_existing <- function(df, cols) {
  hit <- intersect(cols, colnames(df))
  if (length(hit) == 0) {
    stop("Missing required column. Tried: ", paste(cols, collapse = ", "))
  }
  hit[[1]]
}

map_symbols <- function(symbols) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("org.Hs.eg.db and AnnotationDbi are required for Module 7.4 ID mapping")
  }
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  mapped <- mapped[!is.na(mapped$ENTREZID), , drop = FALSE]
  mapped[!duplicated(mapped$SYMBOL), , drop = FALSE]
}

run_ora <- function(tf, genes_entrez, background_entrez, database) {
  if (length(genes_entrez) < 5) {
    return(data.frame())
  }
  if (database == "KEGG") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("clusterProfiler is required for KEGG enrichment")
    }
    res <- clusterProfiler::enrichKEGG(
      gene = genes_entrez,
      universe = background_entrez,
      organism = "hsa",
      pAdjustMethod = "BH"
    )
    table <- as.data.frame(res)
  } else if (database == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("ReactomePA is required for Reactome enrichment")
    }
    res <- ReactomePA::enrichPathway(
      gene = genes_entrez,
      universe = background_entrez,
      organism = "human",
      pAdjustMethod = "BH",
      readable = FALSE
    )
    table <- as.data.frame(res)
  } else {
    stop("Unsupported ORA database: ", database)
  }
  if (nrow(table) == 0) {
    return(data.frame())
  }
  data.frame(
    tf = tf,
    database = database,
    term_id = table$ID,
    term_name = table$Description,
    pvalue = table$pvalue,
    p.adjust = table$p.adjust,
    gene_count = table$Count,
    gene_ratio = table$GeneRatio,
    stringsAsFactors = FALSE
  )
}

run_gsea <- function(tf, ranked_symbols, symbol_to_entrez) {
  if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("fgsea and msigdbr are required for GSEA")
  }
  mapped <- symbol_to_entrez[match(names(ranked_symbols), symbol_to_entrez$SYMBOL), , drop = FALSE]
  keep <- !is.na(mapped$ENTREZID)
  stats <- ranked_symbols[keep]
  names(stats) <- mapped$ENTREZID[keep]
  stats <- sort(tapply(stats, names(stats), max), decreasing = TRUE)
  if (length(stats) < 10) {
    return(data.frame())
  }
  msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  pathways <- split(msig$entrez_gene, msig$gs_name)
  fg <- fgsea::fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500)
  if (nrow(fg) == 0) {
    return(data.frame())
  }
  data.frame(
    tf = tf,
    database = "GSEA_REACTOME",
    term_id = fg$pathway,
    term_name = fg$pathway,
    pvalue = fg$pval,
    p.adjust = fg$padj,
    gene_count = fg$size,
    NES = fg$NES,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args()
  dir.create(args$metadata_dir, recursive = TRUE, showWarnings = FALSE)
  perturb <- read.delim(args$perturbation, stringsAsFactors = FALSE, check.names = FALSE)
  background <- read.delim(args$background_genes, stringsAsFactors = FALSE)[[1]]
  fdr_col <- first_existing(perturb, c("p.adj", "p_adj", "padj", "FDR", "fdr", "qvalue"))
  distance_col <- first_existing(perturb, c("distance", "Distance", "dist", "perturbation_score", "score"))
  perturb$gene <- as.character(perturb$gene)
  perturb[[fdr_col]] <- as.numeric(perturb[[fdr_col]])
  perturb[[distance_col]] <- as.numeric(perturb[[distance_col]])

  symbol_map <- map_symbols(c(background, perturb$gene))
  background_map <- symbol_map[symbol_map$SYMBOL %in% background, , drop = FALSE]
  background_entrez <- unique(background_map$ENTREZID)

  enrichment_rows <- list()
  mapping_rows <- list()
  for (tf in unique(perturb$tf)) {
    tf_table <- perturb[perturb$tf == tf, , drop = FALSE]
    significant <- tf_table[!is.na(tf_table[[fdr_col]]) & tf_table[[fdr_col]] <= args$fdr, , drop = FALSE]
    mapped_input <- symbol_map[symbol_map$SYMBOL %in% significant$gene, , drop = FALSE]
    mapping_rows[[tf]] <- data.frame(
      tf = tf,
      database = "all",
      n_background = length(unique(background)),
      n_input = length(unique(significant$gene)),
      n_mapped_background = length(unique(background_map$SYMBOL)),
      n_mapped_input = length(unique(mapped_input$SYMBOL)),
      stringsAsFactors = FALSE
    )
    genes_entrez <- unique(mapped_input$ENTREZID)
    if (length(genes_entrez) >= args$min_genes) {
      enrichment_rows[[paste(tf, "KEGG")]] <- run_ora(tf, genes_entrez, background_entrez, "KEGG")
      enrichment_rows[[paste(tf, "Reactome")]] <- run_ora(tf, genes_entrez, background_entrez, "Reactome")
    }
    ranked <- -log10(pmax(tf_table[[fdr_col]], .Machine$double.xmin)) * sign(tf_table[[distance_col]])
    names(ranked) <- tf_table$gene
    enrichment_rows[[paste(tf, "GSEA")]] <- run_gsea(tf, ranked, symbol_map)
  }

  enrichment <- dplyr::bind_rows(enrichment_rows)
  mapping <- dplyr::bind_rows(mapping_rows) %>%
    mutate(
      background_mapping_rate = ifelse(n_background > 0, n_mapped_background / n_background, NA_real_),
      input_mapping_rate = ifelse(n_input > 0, n_mapped_input / n_input, NA_real_)
    )

  enrichment_path <- file.path(args$metadata_dir, "sctenifoldknk_module7_4_enrichment_all.tsv")
  mapping_path <- file.path(args$metadata_dir, "sctenifoldknk_module7_4_mapping_stats.tsv")
  report_path <- file.path(args$metadata_dir, "sctenifoldknk_module7_4_r_report.json")
  write.table(enrichment, enrichment_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(mapping, mapping_path, sep = "\t", quote = FALSE, row.names = FALSE)
  jsonlite::write_json(
    list(
      module = "7.4",
      method = "KEGG, Reactome and preranked Reactome GSEA for scTenifoldKnk perturbation genes",
      inputs = list(perturbation = args$perturbation, background_genes = args$background_genes),
      outputs = list(enrichment = enrichment_path, mapping_stats = mapping_path, report = report_path),
      fdr = args$fdr,
      min_genes = args$min_genes
    ),
    report_path,
    pretty = TRUE,
    auto_unbox = TRUE
  )
  message("Wrote Module 7.4 enrichment: ", enrichment_path)
}

main()
