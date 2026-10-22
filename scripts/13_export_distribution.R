#!/usr/bin/env Rscript
# 13_export_distribution.R
# Export per-sample pseudobulk expression distributions for selected genes.
# One row = gene x cell_type x sample. The viewer uses sample-level values for
# violin/box plots; cell-level values are never written.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(optparse)
  library(arrow)
  library(jsonlite)
})

option_list <- list(
  make_option("--input",            type = "character", help = "Seurat .rds"),
  make_option("--output_dir",       type = "character", default = "output/export/distribution"),
  make_option("--prefix",           type = "character", default = "distribution"),
  make_option("--alias",            type = "character", default = ""),
  make_option("--genes",            type = "character", default = ""),
  make_option("--cell_type_column", type = "character", default = "cell_type"),
  make_option("--cell_types",       type = "character", default = ""),
  make_option("--sample_column",    type = "character", default = "sample_id"),
  make_option("--dataset_column",   type = "character", default = "Dataset"),
  make_option("--condition_column", type = "character", default = "condition"),
  make_option("--min_cells",        type = "integer",   default = 10L),
  make_option("--assay",            type = "character", default = "RNA"),
  make_option("--value_layer",      type = "character", default = "data"),
  make_option("--counts_layer",     type = "character", default = "counts"),
  make_option("--exclude_col",      type = "character", default = "Author_Exclude"),
  make_option("--exclude_vals",     type = "character", default = "yes"),
  make_option("--gene_id_type",     type = "character", default = "symbol"),
  make_option("--species",          type = "character", default = "mouse"),
  make_option("--source",           type = "character", default = "HypoMap"),
  make_option("--license",          type = "character", default = "CC BY 4.0"),
  make_option("--citation",         type = "character", default = ""),
  make_option("--notes",            type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(nzchar(opt$input), file.exists(opt$input))
dir.create(file.path(opt$output_dir, "datasets"), recursive = TRUE, showWarnings = FALSE)

split_arg <- function(x) {
  v <- trimws(strsplit(x, ",")[[1]])
  v[nzchar(v)]
}

get_layer <- function(obj, assay, layer) {
  m <- tryCatch(GetAssayData(obj, assay = assay, layer = layer),
                error = function(e) NULL)
  if (is.null(m)) {
    m <- tryCatch(GetAssayData(obj, assay = assay, slot = layer),
                  error = function(e) NULL)
  }
  if (is.null(m)) stop("[distribution] cannot read layer '", layer, "' from assay '", assay, "'")
  m
}

uuid_like <- function(x) {
  f <- tempfile(); writeLines(x, f)
  h <- substr(unname(tools::md5sum(f)), 1, 32); unlink(f)
  paste(substr(h, 1, 8), substr(h, 9, 12), substr(h, 13, 16),
        substr(h, 17, 20), substr(h, 21, 32), sep = "-")
}

message("[distribution] loading ", opt$input)
obj <- readRDS(opt$input)
if (!inherits(obj, "Seurat")) stop("[distribution] input is not a Seurat object")
message("[distribution] ", nrow(obj), " genes x ", ncol(obj), " cells")

md <- obj@meta.data
required <- c(opt$cell_type_column, opt$sample_column)
missing <- setdiff(required, colnames(md))
if (length(missing)) {
  stop("[distribution] missing metadata column(s): ", paste(missing, collapse = ", "))
}
if (!opt$dataset_column %in% colnames(md)) {
  warning("[distribution] dataset_column '", opt$dataset_column, "' absent; dataset will be NA")
}
if (!opt$condition_column %in% colnames(md)) {
  warning("[distribution] condition_column '", opt$condition_column, "' absent; condition will be NA")
}

genes_req <- split_arg(opt$genes)
genes_use <- if (length(genes_req)) {
  missing_genes <- setdiff(genes_req, rownames(obj))
  if (length(missing_genes)) {
    message("[distribution] WARNING requested genes absent: ",
            paste(missing_genes, collapse = ", "))
  }
  intersect(genes_req, rownames(obj))
} else {
  rownames(obj)
}
if (!length(genes_use)) stop("[distribution] no genes to export")

cell_types_req <- split_arg(opt$cell_types)
keep <- !is.na(md[[opt$cell_type_column]]) & !is.na(md[[opt$sample_column]])
if (length(cell_types_req)) {
  keep <- keep & as.character(md[[opt$cell_type_column]]) %in% cell_types_req
}
excl_vals <- split_arg(opt$exclude_vals)
if (nzchar(opt$exclude_col) && opt$exclude_col %in% colnames(md) && length(excl_vals)) {
  keep <- keep & (is.na(md[[opt$exclude_col]]) |
                    !(as.character(md[[opt$exclude_col]]) %in% excl_vals))
}
kept_idx <- which(keep)
if (!length(kept_idx)) stop("[distribution] no cells remain after filtering")

ct_all  <- as.character(md[[opt$cell_type_column]])
sample_all <- as.character(md[[opt$sample_column]])
dataset_all <- if (opt$dataset_column %in% colnames(md))
  as.character(md[[opt$dataset_column]]) else rep(NA_character_, nrow(md))
condition_all <- if (opt$condition_column %in% colnames(md))
  as.character(md[[opt$condition_column]]) else rep(NA_character_, nrow(md))

data_all   <- get_layer(obj, opt$assay, opt$value_layer)
counts_all <- get_layer(obj, opt$assay, opt$counts_layer)

keys <- paste(ct_all[kept_idx], sample_all[kept_idx], sep = " || ")
groups <- split(kept_idx, keys)
rows <- vector("list", length(groups))
n_used <- 0L
n_skipped <- 0L

for (i in seq_along(groups)) {
  idx <- groups[[i]]
  n <- length(idx)
  if (n < opt$min_cells) {
    n_skipped <- n_skipped + 1L
    next
  }
  first <- idx[1]
  d <- data_all[genes_use, idx, drop = FALSE]
  c <- counts_all[genes_use, idx, drop = FALSE]
  mean_expression <- Matrix::rowSums(d) / n
  mean_counts     <- Matrix::rowSums(c) / n
  n_expressing    <- Matrix::rowSums(c > 0)
  pct             <- 100 * n_expressing / n

  rows[[i]] <- data.frame(
    gene_id         = genes_use,
    gene_name       = genes_use,
    cell_type       = ct_all[first],
    sample          = sample_all[first],
    dataset         = dataset_all[first],
    condition       = condition_all[first],
    n_cells         = as.integer(n),
    n_expressing    = as.integer(n_expressing),
    pct_expressing  = as.numeric(pct),
    mean_expression = as.numeric(mean_expression),
    mean_counts     = as.numeric(mean_counts),
    stringsAsFactors = FALSE
  )
  n_used <- n_used + 1L
  if (i %% 50 == 0) gc(verbose = FALSE)
}

out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
if (is.null(out) || !nrow(out)) {
  stop("[distribution] no sample x cell_type groups passed min_cells=", opt$min_cells)
}

gene_index <- do.call(rbind, lapply(genes_use, function(g) {
  sub <- out[out$gene_id == g, , drop = FALSE]
  top <- sub[which.max(sub$mean_expression), , drop = FALSE]
  data.frame(
    gene_id               = g,
    gene_name             = g,
    max_mean_expression   = max(sub$mean_expression, na.rm = TRUE),
    max_pct_expressing    = max(sub$pct_expressing, na.rm = TRUE),
    n_cell_types_detected = length(unique(sub$cell_type[sub$pct_expressing > 0])),
    top_cell_type         = top$cell_type[1],
    top_mean_expression   = top$mean_expression[1],
    stringsAsFactors      = FALSE
  )
}))

ct_summary <- do.call(rbind, lapply(sort(unique(out$cell_type)), function(ct) {
  sub <- out[out$cell_type == ct, , drop = FALSE]
  data.frame(
    cell_type         = ct,
    n_cells           = sum(unique(sub[, c("sample", "n_cells")])$n_cells),
    n_samples         = length(unique(sub$sample)),
    n_genes_detected  = length(unique(sub$gene_id[sub$pct_expressing > 0])),
    n_genes_requested = length(genes_use),
    stringsAsFactors  = FALSE
  )
}))

ds_dir <- file.path(opt$output_dir, "datasets")
pb_file <- paste0(opt$prefix, "_pseudobulk.parquet")
gi_file <- paste0(opt$prefix, "_gene_index.parquet")
ct_file <- paste0(opt$prefix, "_cell_type_summary.parquet")
arrow::write_parquet(out, file.path(ds_dir, pb_file), compression = "zstd")
arrow::write_parquet(gene_index, file.path(ds_dir, gi_file), compression = "zstd")
arrow::write_parquet(ct_summary, file.path(ds_dir, ct_file), compression = "zstd")

organism <- if (opt$species == "mouse") "Mus musculus" else "Homo sapiens"
alias <- if (nzchar(opt$alias)) opt$alias else paste(opt$source, "pseudobulk distribution")
notes <- if (nzchar(opt$notes)) opt$notes else
  paste0(opt$source, "; per-sample pseudobulk; ", length(genes_use), " selected genes")

item <- list(
  dataset_id        = uuid_like(pb_file),
  alias             = alias,
  original_filename = pb_file,
  file_path         = pb_file,
  dataset_type      = "expression_distribution",
  organism          = organism,
  row_count         = nrow(out),
  gene_count        = length(genes_use),
  cell_type_count   = length(unique(out$cell_type)),
  sample_count      = length(unique(out$sample)),
  notes             = notes,
  tags              = as.list(c("scRNA", "expression", "distribution", "pseudobulk")),
  distribution_schema = list(
    schema_version   = "1.0",
    layout           = "pseudobulk",
    gene_id_type     = opt$gene_id_type,
    cell_type_column = "cell_type",
    sample_column    = "sample",
    sample_unit      = "biological_replicate",
    dataset_column   = "dataset",
    condition_column = "condition",
    group_columns    = as.list(c("condition", "dataset")),
    value_column     = "mean_expression",
    pct_column       = "pct_expressing",
    n_cells_column   = "n_cells",
    pct_scale        = "0-100",
    detection        = "raw counts > 0",
    value_definition = "mean of normalized value layer per sample x cell_type",
    n_genes          = length(genes_use),
    n_samples        = length(unique(out$sample)),
    min_cells        = opt$min_cells,
    source           = opt$source,
    license          = opt$license,
    citation         = opt$citation,
    companion_files  = list(gene_index = gi_file, cell_type_summary = ct_file)
  )
)

jsonlite::write_json(
  list(schema_version = "1.0", datasets = list(item)),
  file.path(opt$output_dir, "seqviewer_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, null = "null"
)

columns_md <- c(
  "# expression_distribution (per-sample pseudobulk)",
  "",
  paste0("Source: `", normalizePath(opt$input), "`"),
  paste0("Layout: `pseudobulk` | sample unit: `biological_replicate` | min_cells: ", opt$min_cells),
  "",
  "One row = gene x cell_type x sample. Values are aggregated; no cell-level data.",
  "",
  "| column | description |",
  "|---|---|",
  "| gene_id, gene_name | gene symbol |",
  "| cell_type | cell type label |",
  "| sample | biological sample/replicate ID |",
  "| dataset | source dataset/study |",
  "| condition | experimental condition |",
  "| n_cells | cells in this sample x cell_type group |",
  "| n_expressing | cells with raw counts > 0 |",
  "| pct_expressing | 100 * n_expressing / n_cells |",
  "| mean_expression | mean normalized value layer |",
  "| mean_counts | mean raw counts |",
  "",
  "Companion files:",
  paste0("- `datasets/", gi_file, "`: per-gene summary"),
  paste0("- `datasets/", ct_file, "`: per-cell-type summary"),
  ""
)
writeLines(columns_md, file.path(opt$output_dir, "COLUMNS.md"))

message("[distribution] ", nrow(out), " rows | ", length(genes_use), " genes | ",
        length(unique(out$cell_type)), " cell types | ", length(unique(out$sample)),
        " samples | used ", n_used, " groups, skipped ", n_skipped,
        " -> ", opt$output_dir)
