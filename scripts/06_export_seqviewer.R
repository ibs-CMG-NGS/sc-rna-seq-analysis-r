#!/usr/bin/env Rscript
# 06_export_seqviewer.R
# Export pipeline results as a cmg-seqviewer handoff bundle:
# parquet datasets + seqviewer_manifest.json + export_index.tsv + COLUMNS.md.
#
# The pipeline is the producer; cmg-seqviewer is the consumer. The viewer never
# reads Seurat/.h5ad objects. Column names and manifest keys follow the viewer
# team's reply (docs/developer/pipeline-integration-response.md).

suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option("--de_dir",          type = "character", help = "Directory with de_all_celltypes.csv and pseudobulk_counts.rds"),
  make_option("--enrichment_dir",  type = "character", default = "", help = "Enrichment CSV directory (optional)"),
  make_option("--composition_dir", type = "character", default = "", help = "Composition directory with propeller_results.csv (optional)"),
  make_option("--output_dir",      type = "character", default = "output/export/seqviewer", help = "Export root [default: %default]"),
  make_option("--schema_version",  type = "character", default = "1.0", help = "Export schema version [default: %default]"),
  make_option("--species",         type = "character", default = "human", help = "human | mouse [default: %default]"),
  make_option("--padj_cutoff",     type = "double",    default = 0.05, help = "[default: %default]"),
  make_option("--log2fc_cutoff",   type = "double",    default = 1.0,  help = "[default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$de_dir)) stop("--de_dir is required", call. = FALSE)

for (pkg in c("arrow", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("[export_seqviewer] R package '", pkg, "' is required. ",
         "Install with: conda install -c conda-forge r-", pkg, call. = FALSE)
  }
}

organism <- if (opt$species == "mouse") "Mus musculus" else "Homo sapiens"
safe     <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

datasets_dir <- file.path(opt$output_dir, "datasets")
dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(datasets_dir, pattern = "\\.parquet$", full.names = TRUE))

entries    <- list()
index_rows <- list()

uuid5 <- function(name, ns = "6ba7b810-9dad-11d1-80b4-00c04fd430c8") {
  ns_bytes <- as.raw(strtoi(substring(gsub("-", "", ns), seq(1, 31, 2), seq(2, 32, 2)), 16L))
  h  <- digest::digest(c(ns_bytes, charToRaw(name)), algo = "sha1", serialize = FALSE, raw = TRUE)
  b  <- as.integer(h[1:16])
  b[7] <- bitwOr(bitwAnd(b[7], 0x0FL), 0x50L)
  b[9] <- bitwOr(bitwAnd(b[9], 0x3FL), 0x80L)
  hx <- sprintf("%02x", b)
  paste0(paste(hx[1:4], collapse = ""), "-", paste(hx[5:6], collapse = ""), "-",
         paste(hx[7:8], collapse = ""), "-", paste(hx[9:10], collapse = ""), "-",
         paste(hx[11:16], collapse = ""))
}

col_or <- function(df, name, type = "chr") {
  if (name %in% names(df)) {
    return(if (type == "num") as.numeric(df[[name]]) else as.character(df[[name]]))
  }
  if (type == "num") rep(NA_real_, nrow(df)) else rep(NA_character_, nrow(df))
}

register <- function(df, stem, dataset_type, alias, condition, cell_type = NA_character_,
                     direction = NA_character_, ontology = NA_character_,
                     gene_count = 0L, significant = 0L, notes = "", tags = character(0)) {
  fname <- paste0(stem, ".parquet")
  arrow::write_parquet(df, file.path(datasets_dir, fname))

  entry <- list(
    dataset_id           = uuid5(fname),
    alias                = alias,
    original_filename    = fname,
    file_path            = fname,
    dataset_type         = dataset_type,
    experiment_condition = condition,
    cell_type            = cell_type,
    organism             = organism,
    row_count            = as.integer(nrow(df)),
    gene_count           = as.integer(gene_count),
    significant_genes    = as.integer(significant),
    analysis_date        = format(Sys.Date()),
    notes                = notes,
    tags                 = as.list(unique(c("scRNA", tags)))
  )
  is_absent <- vapply(entry, function(x) length(x) == 1 && !is.list(x) && is.na(x), logical(1))
  entries[[length(entries) + 1]] <<- entry[!is_absent]

  index_rows[[length(index_rows) + 1]] <<- data.frame(
    file = fname, dataset_type = dataset_type, cell_type = cell_type,
    condition = condition, direction = direction, ontology = ontology,
    rows = nrow(df), stringsAsFactors = FALSE
  )
  message("[export_seqviewer] ", sprintf("%-18s", dataset_type), " ", fname, ": ", nrow(df), " rows")
}

# ── DE: one dataset per cell_type x contrast ─────────────────────────────────
de_file <- file.path(opt$de_dir, "de_all_celltypes.csv")
if (!file.exists(de_file)) stop("[export_seqviewer] missing ", de_file, call. = FALSE)
de <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)

need <- c("gene", "baseMean", "log2FoldChange", "pvalue", "padj", "contrast", "cell_type")
miss <- setdiff(need, names(de))
if (length(miss) > 0) stop("[export_seqviewer] DE table missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

de_std <- data.frame(
  gene_id    = as.character(de$gene),
  symbol     = as.character(de$gene),
  base_mean  = as.numeric(de$baseMean),
  log2fc     = as.numeric(de$log2FoldChange),
  lfcse      = col_or(de, "lfcSE", "num"),
  stat       = col_or(de, "stat", "num"),
  pvalue     = as.numeric(de$pvalue),
  adj_pvalue = as.numeric(de$padj),
  cell_type  = as.character(de$cell_type),
  contrast   = as.character(de$contrast),
  stringsAsFactors = FALSE
)

de_pairs <- unique(de_std[, c("cell_type", "contrast")])
for (i in seq_len(nrow(de_pairs))) {
  ct  <- de_pairs$cell_type[i]
  ctr <- de_pairs$contrast[i]
  sub <- de_std[de_std$cell_type == ct & de_std$contrast == ctr, , drop = FALSE]
  sig <- sum(sub$adj_pvalue < opt$padj_cutoff & abs(sub$log2fc) >= opt$log2fc_cutoff, na.rm = TRUE)
  register(
    sub, stem = paste0("DE__", safe(ct), "__", safe(ctr)),
    dataset_type = "differential_expression",
    alias = paste(ct, ctr, "DE"), condition = ctr, cell_type = ct,
    gene_count = nrow(sub), significant = sig,
    notes = sprintf("DESeq2 pseudobulk; padj < %g, |log2FC| >= %g; schema %s",
                    opt$padj_cutoff, opt$log2fc_cutoff, opt$schema_version),
    tags = c("DE", "DESeq2", ct, ctr)
  )
}

# ── GO / KEGG: one dataset per cell_type x contrast x database x direction ───
if (nzchar(opt$enrichment_dir) && dir.exists(opt$enrichment_dir)) {
  label_keys <- data.frame(
    label_safe = safe(paste0(de_pairs$cell_type, "_", de_pairs$contrast)),
    cell_type  = de_pairs$cell_type,
    contrast   = de_pairs$contrast,
    stringsAsFactors = FALSE
  )
  label_keys <- label_keys[order(-nchar(label_keys$label_safe)), ]

  enr_files <- list.files(opt$enrichment_dir, pattern = "\\.csv$", full.names = TRUE)
  for (f in enr_files) {
    base <- sub("\\.csv$", "", basename(f))
    key  <- Find(function(k) startsWith(base, paste0(k, "_")), label_keys$label_safe)
    if (is.null(key)) {
      warning("[export_seqviewer] could not match enrichment file: ", basename(f))
      next
    }
    key_row <- label_keys[label_keys$label_safe == key, ]
    db_name <- sub(paste0("^", key, "_"), "", base)
    tok     <- strsplit(db_name, "_")[[1]]

    if (identical(tok[1], "GO") && length(tok) == 3) {
      ontology <- tok[2]; direction <- toupper(tok[3])
    } else if (identical(tok[1], "KEGG") && length(tok) == 2) {
      ontology <- "KEGG"; direction <- toupper(tok[2])
    } else {
      warning("[export_seqviewer] unrecognized enrichment suffix: ", db_name)
      next
    }

    df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0) next
    if (!all(c("ID", "Description", "Count", "p.adjust") %in% names(df))) {
      warning("[export_seqviewer] enrichment table missing columns: ", basename(f))
      next
    }

    out <- data.frame(
      term_id      = as.character(df$ID),
      description  = as.character(df$Description),
      gene_count   = as.integer(df$Count),
      gene_ratio   = col_or(df, "GeneRatio"),
      bg_ratio     = col_or(df, "BgRatio"),
      pvalue       = col_or(df, "pvalue", "num"),
      fdr          = as.numeric(df$p.adjust),
      qvalue       = col_or(df, "qvalue", "num"),
      gene_symbols = col_or(df, "geneID"),
      gene_set     = direction,
      ontology     = ontology,
      direction    = direction,
      cell_type    = key_row$cell_type,
      contrast     = key_row$contrast,
      stringsAsFactors = FALSE
    )

    register(
      out, stem = paste0("ENR__", key, "__", db_name),
      dataset_type = "go_analysis",
      alias = paste(key_row$cell_type, key_row$contrast, db_name),
      condition = key_row$contrast, cell_type = key_row$cell_type,
      direction = direction, ontology = ontology,
      gene_count = 0L, significant = 0L,
      notes = sprintf("clusterProfiler ORA; q < 0.2 pipeline cutoff; schema %s", opt$schema_version),
      tags = c("enrichment", ontology, direction, key_row$cell_type, key_row$contrast)
    )
  }
}

# ── Composition: propeller result table ───────────────────────────────────────
if (nzchar(opt$composition_dir)) {
  prop_file <- file.path(opt$composition_dir, "propeller_results.csv")
  if (file.exists(prop_file)) {
    pr <- read.csv(prop_file, check.names = FALSE, stringsAsFactors = FALSE)
    pm <- grep("^PropMean\\.", names(pr), value = TRUE)
    keep <- c("cell_type", pm, intersect(c("PropRatio", "Tstatistic", "P.Value", "FDR"), names(pr)))
    out <- pr[, keep, drop = FALSE]
    conds <- sub("^PropMean\\.", "", pm)
    register(
      out, stem = "COMP__propeller",
      dataset_type = "composition",
      alias = "Cell-type composition (propeller)",
      condition = paste(conds, collapse = " vs "),
      gene_count = nrow(out),
      significant = sum(out$FDR < opt$padj_cutoff, na.rm = TRUE),
      notes = sprintf("propeller (speckle); FDR < %g; schema %s", opt$padj_cutoff, opt$schema_version),
      tags = c("composition", "propeller")
    )
  } else {
    message("[export_seqviewer] no propeller_results.csv; composition skipped")
  }
}

# ── Multi-Group: pseudobulk counts per cell type (samples as <condition>_<rep>) ─
pb_file <- file.path(opt$de_dir, "pseudobulk_counts.rds")
if (file.exists(pb_file)) {
  pb <- readRDS(pb_file)
  for (ct in names(pb$col_data)) {
    cd   <- pb$col_data[[ct]]
    srcs <- paste0(ct, "__", rownames(cd))
    if (!all(srcs %in% colnames(pb$counts)) || nrow(cd) < 3) {
      warning("[export_seqviewer] skipping Multi-Group for ", ct, ": insufficient samples")
      next
    }
    cd$rep    <- ave(seq_len(nrow(cd)), cd$condition, FUN = seq_along)
    mat       <- as.matrix(pb$counts[, srcs, drop = FALSE])
    out       <- data.frame(gene_id = rownames(pb$counts), mat,
                            check.names = FALSE, stringsAsFactors = FALSE)
    colnames(out)[-1] <- paste0(cd$condition, "_", cd$rep)

    register(
      out, stem = paste0("MG__", safe(ct)),
      dataset_type = "multi_group",
      alias = paste(ct, "pseudobulk counts"),
      condition = paste(unique(cd$condition), collapse = " vs "),
      cell_type = ct,
      gene_count = nrow(out),
      notes = sprintf("pseudobulk raw counts (UMI); schema %s", opt$schema_version),
      tags = c("pseudobulk", "counts", ct)
    )
  }
}

# ── Manifest, index, column contract ─────────────────────────────────────────
jsonlite::write_json(list(datasets = entries),
                     file.path(opt$output_dir, "seqviewer_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE)

index <- do.call(rbind, index_rows)
index <- index[order(index$dataset_type, index$cell_type, index$condition), ]
write.table(index, file.path(opt$output_dir, "export_index.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

columns_md <- c(
  "# cmg-seqviewer handoff: column contract",
  "",
  paste0("schema_version: ", opt$schema_version),
  "",
  "Producer: scRNA-seq 3rd analysis pipeline. Consumer: cmg-seqviewer.",
  "The bundle root is an Import Folder: `seqviewer_manifest.json` + `datasets/*.parquet`.",
  "",
  "## differential_expression (one file per cell_type x contrast)",
  "",
  "| column | required | source | notes |",
  "| --- | --- | --- | --- |",
  "| gene_id | yes | gene | gene symbol (no Ensembl ID in current input) |",
  "| symbol | yes | gene | same as gene_id |",
  "| log2fc | yes | log2FoldChange | |",
  "| adj_pvalue | yes | padj | BH-adjusted |",
  "| base_mean | no | baseMean | |",
  "| lfcse | no | lfcSE | |",
  "| stat | no | stat | |",
  "| pvalue | no | pvalue | |",
  "| cell_type | no | cell_type | metadata |",
  "| contrast | no | contrast | metadata |",
  "",
  "## go_analysis (one file per cell_type x contrast x database x direction)",
  "",
  "| column | required | source | notes |",
  "| --- | --- | --- | --- |",
  "| term_id | yes | ID | GO or KEGG id |",
  "| description | yes | Description | |",
  "| gene_count | yes | Count | |",
  "| fdr | yes | p.adjust | BH-adjusted (not `padj`) |",
  "| gene_ratio, bg_ratio | no | GeneRatio, BgRatio | |",
  "| pvalue, qvalue | no | pvalue, qvalue | |",
  "| gene_symbols | no | geneID | `/`-separated |",
  "| gene_set | no | direction | UP / DOWN |",
  "| ontology | no | derived | BP / CC / MF / KEGG |",
  "| direction | no | derived | UP / DOWN |",
  "| cell_type, contrast | no | derived | metadata |",
  "",
  "## composition (one file)",
  "",
  "| column | required | notes |",
  "| --- | --- | --- |",
  "| cell_type | yes | |",
  "| PropMean.<condition> | yes (or PropRatio) | per-condition mean proportion |",
  "| PropRatio | yes (or 2+ PropMean) | |",
  "| P.Value, FDR | recommended | propeller output |",
  "",
  "## multi_group (one file per cell type)",
  "",
  "| column | notes |",
  "| --- | --- |",
  "| gene_id | first column |",
  "| <condition>_<rep> | raw counts, numeric; at least 3 samples |",
  "",
  "`seqviewer_manifest.json` sets `dataset_type` explicitly for every entry.",
  "`export_index.tsv` lists the files for pipeline-side checks."
)
writeLines(columns_md, file.path(opt$output_dir, "COLUMNS.md"))

message("[export_seqviewer] Done: ", length(entries), " dataset(s) -> ", opt$output_dir)
