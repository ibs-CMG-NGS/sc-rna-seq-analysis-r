#!/usr/bin/env Rscript
# 06_export_seqviewer.R
# Export pseudobulk DE and enrichment results as cmg-seqviewer-ready parquet
# handoff bundles.
#
# Design contract (v3.0): the pipeline is the *producer*; cmg-seqviewer is the
# *consumer*. The viewer never reads Seurat/.h5ad objects directly — it consumes
# plain tables reshaped into its existing DE / GO·KEGG schemas. Output is
# tool-agnostic (parquet + a column contract) so any table-based viewer can reuse
# it. See docs/cmg-seqviewer-integration-questions.md and .claude/PLAN.md (Track B).

suppressPackageStartupMessages({
  library(optparse)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--de_dir",         type = "character",
              help = "Directory containing de_*.csv files"),
  make_option("--enrichment_dir", type = "character", default = "",
              help = "Directory containing enrichment *.csv files (optional)"),
  make_option("--output_dir",     type = "character", default = "output/export/seqviewer",
              help = "Export root [default: %default]"),
  make_option("--schema_version", type = "character", default = "1.0",
              help = "Column-contract version recorded in the index [default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$de_dir))

if (!requireNamespace("arrow", quietly = TRUE)) {
  stop("[export_seqviewer] R package 'arrow' is required for parquet export.\n",
       "  Install with: conda install -c conda-forge r-arrow", call. = FALSE)
}

datasets_dir <- file.path(opt$output_dir, "datasets")
dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)

index_rows <- list()
add_index <- function(file, dataset_type, cell_type = NA_character_,
                      contrast = NA_character_, direction = NA_character_,
                      ontology = NA_character_, rows = NA_integer_) {
  index_rows[[length(index_rows) + 1]] <<- data.frame(
    file = file, dataset_type = dataset_type, cell_type = cell_type,
    contrast = contrast, direction = direction, ontology = ontology,
    rows = rows, stringsAsFactors = FALSE
  )
}

write_parquet <- function(df, path) {
  arrow::write_parquet(df, path)
  path
}

# ── Column mapping ────────────────────────────────────────────────────────────
# cmg-seqviewer canonical names on the left; our internal names on the right.
# Extra descriptive columns (lfcSE/stat/pvalue/cell_type/contrast) are kept so a
# single export stays self-describing (see Q3 in the integration questions doc).
map_de <- function(df) {
  need <- c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) {
    warning("[export_seqviewer] DE table missing columns: ", paste(miss, collapse = ", "))
    return(NULL)
  }
  data.frame(
    gene_id    = as.character(df$gene),
    base_mean  = as.numeric(df$baseMean),
    log2fc     = as.numeric(df$log2FoldChange),
    lfcSE      = if ("lfcSE" %in% names(df)) as.numeric(df$lfcSE) else NA_real_,
    stat       = if ("stat"  %in% names(df)) as.numeric(df$stat)  else NA_real_,
    pvalue     = as.numeric(df$pvalue),
    adj_pvalue = as.numeric(df$padj),
    cell_type  = if ("cell_type" %in% names(df)) as.character(df$cell_type) else NA_character_,
    contrast   = if ("contrast"  %in% names(df)) as.character(df$contrast)  else NA_character_,
    stringsAsFactors = FALSE
  )
}

map_enrichment <- function(df, db, ontology, direction, cell_type, contrast) {
  need <- c("ID", "Description", "Count", "p.adjust")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) {
    warning("[export_seqviewer] enrichment table missing columns: ", paste(miss, collapse = ", "))
    return(NULL)
  }
  n        <- nrow(df)
  gene_set <- if (identical(db, "GO")) paste0("GO:", ontology) else "KEGG"
  data.frame(
    term_id         = as.character(df$ID),
    description     = as.character(df$Description),
    gene_count      = as.integer(df$Count),
    gene_ratio      = if ("GeneRatio"      %in% names(df)) as.character(df$GeneRatio) else rep(NA_character_, n),
    fold_enrichment = if ("FoldEnrichment" %in% names(df)) as.numeric(df$FoldEnrichment) else rep(NA_real_, n),
    pvalue          = if ("pvalue"         %in% names(df)) as.numeric(df$pvalue) else rep(NA_real_, n),
    fdr             = as.numeric(df$p.adjust),
    qvalue          = if ("qvalue"         %in% names(df)) as.numeric(df$qvalue) else rep(NA_real_, n),
    gene_symbols    = if ("geneID"         %in% names(df)) as.character(df$geneID) else rep(NA_character_, n),
    gene_set        = rep(gene_set, n),
    ontology        = rep(ontology, n),
    direction       = rep(direction, n),
    cell_type       = rep(cell_type, n),
    contrast        = rep(contrast, n),
    stringsAsFactors = FALSE
  )
}

# ── DE export ─────────────────────────────────────────────────────────────────
de_files <- list.files(opt$de_dir, pattern = "^de_.*\\.csv$", full.names = TRUE)
per_type <- de_files[!grepl("de_(all_celltypes|significant)\\.csv$", de_files)]
combined <- de_files[grepl("de_all_celltypes\\.csv$", de_files)]

de_contrasts <- list()   # cell_type -> contrasts (for enrichment key reconstruction)
for (f in per_type) {
  ct <- sub("^de_(.+)\\.csv$", "\\1", basename(f))
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  out <- map_de(df)
  if (is.null(out)) next
  de_contrasts[[ct]] <- unique(stats::na.omit(out$contrast))

  path <- file.path(datasets_dir, paste0("DE__", ct, ".parquet"))
  write_parquet(out, path)
  add_index(basename(path), "DE", cell_type = ct,
            contrast = paste(de_contrasts[[ct]], collapse = ";"), rows = nrow(out))
  message("[export_seqviewer] DE  ", ct, ": ", nrow(out), " rows -> ", basename(path))
}

if (length(combined) > 0) {
  df  <- read.csv(combined[1], stringsAsFactors = FALSE, check.names = FALSE)
  out <- map_de(df)
  if (!is.null(out)) {
    path <- file.path(datasets_dir, "DE__all.parquet")
    write_parquet(out, path)
    add_index(basename(path), "DE", cell_type = "all",
              contrast = paste(unique(out$contrast), collapse = ";"), rows = nrow(out))
    message("[export_seqviewer] DE  all: ", nrow(out), " rows -> ", basename(path))
  }
}

# ── Enrichment export ─────────────────────────────────────────────────────────
if (nzchar(opt$enrichment_dir) && dir.exists(opt$enrichment_dir)) {
  enr_files <- list.files(opt$enrichment_dir, pattern = "\\.csv$", full.names = TRUE)
  # Reconstruct {cell_type}_{contrast} label_safe keys from the DE side so the
  # enrichment filename (which concatenates them, both underscore-separated) is
  # never parsed ambiguously.
  key_map <- list()
  for (ct in names(de_contrasts)) {
    for (ctr in de_contrasts[[ct]]) {
      label_safe <- gsub("[^A-Za-z0-9_]", "_", paste0(ct, "_", ctr))
      key_map[[label_safe]] <- list(cell_type = ct, contrast = ctr)
    }
  }

  for (f in enr_files) {
    base <- sub("\\.csv$", "", basename(f))
    hit <- NULL
    for (ks in names(key_map)) {
      if (startsWith(base, paste0(ks, "_"))) { hit <- ks; break }
    }
    if (is.null(hit)) {
      warning("[export_seqviewer] could not match enrichment file: ", basename(f))
      next
    }
    db_name <- sub(paste0("^", hit, "_"), "", base)   # e.g. GO_BP_up | KEGG_up
    tok     <- strsplit(db_name, "_")[[1]]
    db      <- tok[1]
    if (identical(db, "GO") && length(tok) == 3) {
      ontology <- tok[2]; direction <- tok[3]
    } else if (identical(db, "KEGG") && length(tok) == 2) {
      ontology <- "KEGG"; direction <- tok[2]
    } else {
      warning("[export_seqviewer] unrecognized enrichment db suffix: ", db_name)
      next
    }

    df  <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(df) == 0) {
      message("[export_seqviewer] skip empty enrichment table: ", basename(f))
      next
    }
    out <- map_enrichment(df, db, ontology, direction,
                          key_map[[hit]]$cell_type, key_map[[hit]]$contrast)
    if (is.null(out)) next

    path <- file.path(datasets_dir, paste0("ENR__", hit, "__", db_name, ".parquet"))
    write_parquet(out, path)
    add_index(basename(path), db, cell_type = key_map[[hit]]$cell_type,
              contrast = key_map[[hit]]$contrast, direction = direction,
              ontology = ontology, rows = nrow(out))
    message("[export_seqviewer] ", sprintf("%-4s", db), hit, " ", db_name,
            ": ", nrow(out), " rows -> ", basename(path))
  }
}

# ── Index + column contract ───────────────────────────────────────────────────
index <- if (length(index_rows) > 0) do.call(rbind, index_rows) else
  data.frame(file = character(), dataset_type = character(), cell_type = character(),
             contrast = character(), direction = character(), ontology = character(),
             rows = integer(), stringsAsFactors = FALSE)

index <- index[order(index$dataset_type, index$cell_type, index$contrast, index$direction), ]
idx_path <- file.path(opt$output_dir, "export_index.tsv")
write.table(index, idx_path, sep = "\t", row.names = FALSE, quote = FALSE)

columns_md <- c(
  "# cmg-seqviewer handoff — column contract",
  "",
  paste0("schema_version: ", opt$schema_version),
  "",
  "> Producer: scRNA-seq 3rd analysis pipeline (R/Seurat).",
  "> Consumer: cmg-seqviewer. The viewer never reads Seurat/.h5ad objects; it",
  "> consumes the parquet tables below, reshaped into its existing schemas.",
  "> Tool-agnostic: any table-based viewer can reuse these files.",
  "",
  "## DE dataset (dataset_type = DE)",
  "",
  "| column | viewer-canonical | source (internal) | notes |",
  "| --- | --- | --- | --- |",
  "| gene_id | yes | gene | gene symbol |",
  "| log2fc | yes | log2FoldChange | |",
  "| adj_pvalue | yes | padj | BH-adjusted |",
  "| base_mean | yes | baseMean | mean normalized count |",
  "| pvalue | no | pvalue | kept for context |",
  "| lfcSE, stat | no | lfcSE, stat | kept for context |",
  "| cell_type | no | cell_type | metadata (Q3) |",
  "| contrast | no | contrast | metadata (Q3), e.g. condition_treat_vs_ctrl |",
  "",
  "## GO / KEGG dataset (dataset_type = GO | KEGG)",
  "",
  "| column | viewer-canonical | source (internal) | notes |",
  "| --- | --- | --- | --- |",
  "| term_id | yes | ID | GO id / KEGG pathway id |",
  "| description | yes | Description | |",
  "| gene_count | yes | Count | |",
  "| fdr | yes | p.adjust | BH-adjusted (NOT 'padj') |",
  "| gene_ratio | yes | GeneRatio | e.g. 5/200 |",
  "| fold_enrichment | yes | FoldEnrichment | |",
  "| gene_symbols | yes | geneID | '/' separated symbols |",
  "| gene_set | yes | derived | GO:BP/GO:CC/GO:MF or KEGG |",
  "| ontology | yes | derived | BP/CC/MF or KEGG (Q5: confirm expected encoding) |",
  "| pvalue, qvalue | no | pvalue, qvalue | kept for context |",
  "| direction | no | derived | up / down |",
  "| cell_type, contrast | no | derived | metadata |",
  "",
  "Datasets live in `datasets/*.parquet`; `export_index.tsv` lists them.",
  "`seqviewer_manifest.json` (DB import metadata) is produced by a later phase (E2)."
)
writeLines(columns_md, file.path(opt$output_dir, "COLUMNS.md"))

message("[export_seqviewer] Done. ", nrow(index), " dataset(s) in ", opt$output_dir)
message("[export_seqviewer]   index : ", idx_path)
message("[export_seqviewer]   schema: ", file.path(opt$output_dir, "COLUMNS.md"))
