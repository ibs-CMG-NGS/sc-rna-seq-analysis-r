#!/usr/bin/env Rscript
# 02_enrichment.R
# GO and KEGG enrichment on pseudobulk DE results (per cell type, per contrast).

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--de_dir",         type = "character", help = "Directory containing de_*.csv files"),
  make_option("--species",        type = "character", default = "human",
              help = "Species: human | mouse [default: %default]"),
  make_option("--go_ontologies",  type = "character", default = "BP,CC,MF",
              help = "Comma-separated GO ontologies [default: %default]"),
  make_option("--run_kegg",       type = "character", default = "TRUE",
              help = "Run KEGG enrichment: TRUE | FALSE [default: %default]"),
  make_option("--pvalue_cutoff",  type = "double",    default = 0.05),
  make_option("--qvalue_cutoff",  type = "double",    default = 0.2),
  make_option("--min_gs_size",    type = "integer",   default = 10L),
  make_option("--max_gs_size",    type = "integer",   default = 500L),
  make_option("--output_dir",     type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$de_dir), !is.null(opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
run_kegg    <- toupper(opt$run_kegg) == "TRUE"
go_ontos    <- trimws(strsplit(opt$go_ontologies, ",")[[1]])

# ── OrgDb selection ───────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
})

if (opt$species == "human") {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  orgdb    <- org.Hs.eg.db
  kegg_org <- "hsa"
} else if (opt$species == "mouse") {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  orgdb    <- org.Mm.eg.db
  kegg_org <- "mmu"
} else {
  stop("Unknown species: ", opt$species, ". Use 'human' or 'mouse'.")
}

# ── Load DE results ───────────────────────────────────────────────────────────
de_files <- list.files(opt$de_dir, pattern = "^de_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("de_all_celltypes|de_significant", de_files)]

if (length(de_files) == 0) {
  # Fall back to combined file
  combined_file <- file.path(opt$de_dir, "de_all_celltypes.csv")
  if (!file.exists(combined_file)) stop("No DE result files found in: ", opt$de_dir)
  combined <- read.csv(combined_file)
  de_list  <- split(combined, combined$cell_type)
} else {
  de_list <- lapply(de_files, read.csv)
  names(de_list) <- sub("^de_(.+)\\.csv$", "\\1", basename(de_files))
}

message("[enrichment] Processing ", length(de_list), " cell type(s)...")

# ── Enrichment helper ─────────────────────────────────────────────────────────
run_enrichment <- function(genes_up, genes_dn, all_genes, label, plot_dir) {

  results <- list()

  for (direction in c("up", "down")) {
    gene_set <- if (direction == "up") genes_up else genes_dn
    if (length(gene_set) < 3) next

    # Gene name -> Entrez
    eg <- tryCatch(
      bitr(gene_set, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb),
      error = function(e) NULL
    )
    if (is.null(eg) || nrow(eg) == 0) next
    entrez_genes <- eg$ENTREZID

    # Universe
    universe_eg <- tryCatch(
      bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)$ENTREZID,
      error = function(e) NULL
    )

    # GO
    for (onto in go_ontos) {
      go_res <- tryCatch(
        enrichGO(
          gene          = entrez_genes,
          universe      = universe_eg,
          OrgDb         = orgdb,
          ont           = onto,
          pAdjustMethod = "BH",
          pvalueCutoff  = opt$pvalue_cutoff,
          qvalueCutoff  = opt$qvalue_cutoff,
          minGSSize     = opt$min_gs_size,
          maxGSSize     = opt$max_gs_size,
          readable      = TRUE
        ),
        error = function(e) NULL
      )
      if (!is.null(go_res) && nrow(go_res@result) > 0) {
        key <- paste0("GO_", onto, "_", direction)
        results[[key]] <- as.data.frame(go_res)

        # Dotplot
        p <- dotplot(go_res, showCategory = 20) +
          labs(title = paste0(label, "\n", direction, " - GO:", onto)) +
          theme(plot.title = element_text(size = 9))
        label_safe <- gsub("[^A-Za-z0-9_]", "_", label)
        ggsave(file.path(plot_dir, paste0(label_safe, "_GO_", onto, "_", direction, ".png")),
               p, width = 8, height = 7, dpi = 150)
      }
    }

    # KEGG
    if (run_kegg) {
      kegg_res <- tryCatch(
        enrichKEGG(
          gene          = entrez_genes,
          universe      = universe_eg,
          organism      = kegg_org,
          pAdjustMethod = "BH",
          pvalueCutoff  = opt$pvalue_cutoff,
          qvalueCutoff  = opt$qvalue_cutoff,
          minGSSize     = opt$min_gs_size,
          maxGSSize     = opt$max_gs_size
        ),
        error = function(e) NULL
      )
      if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
        key <- paste0("KEGG_", direction)
        results[[key]] <- as.data.frame(kegg_res)

        p <- dotplot(kegg_res, showCategory = 20) +
          labs(title = paste0(label, "\n", direction, " - KEGG")) +
          theme(plot.title = element_text(size = 9))
        label_safe <- gsub("[^A-Za-z0-9_]", "_", label)
        ggsave(file.path(plot_dir, paste0(label_safe, "_KEGG_", direction, ".png")),
               p, width = 8, height = 7, dpi = 150)
      }
    }
  }
  results
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_enr <- list()

for (ct_name in names(de_list)) {
  df <- de_list[[ct_name]]
  if (!all(c("gene", "padj", "log2FoldChange") %in% colnames(df))) {
    warning("Skipping ", ct_name, ": missing required columns.")
    next
  }

  contrasts <- unique(df$contrast)
  ct_safe   <- gsub("[^A-Za-z0-9_]", "_", ct_name)
  ct_plot_dir <- file.path(opt$output_dir, "plots", ct_safe)
  dir.create(ct_plot_dir, recursive = TRUE, showWarnings = FALSE)

  for (ctr in contrasts) {
    sub     <- filter(df, contrast == ctr, !is.na(padj), !is.na(log2FoldChange))
    genes_up <- sub$gene[sub$padj < opt$pvalue_cutoff & sub$log2FoldChange > 0]
    genes_dn <- sub$gene[sub$padj < opt$pvalue_cutoff & sub$log2FoldChange < 0]
    all_genes <- sub$gene

    label  <- paste0(ct_name, "_", ctr)
    message("  ", label, ": up=", length(genes_up), ", down=", length(genes_dn))

    enr <- run_enrichment(genes_up, genes_dn, all_genes, label, ct_plot_dir)
    all_enr[[label]] <- enr

    # Save CSV per cell_type × contrast × database
    for (db_name in names(enr)) {
      label_safe <- gsub("[^A-Za-z0-9_]", "_", label)
      write.csv(enr[[db_name]],
                file = file.path(opt$output_dir, paste0(label_safe, "_", db_name, ".csv")),
                row.names = FALSE)
    }
  }
}

message("[enrichment] Done. Results in: ", opt$output_dir)
