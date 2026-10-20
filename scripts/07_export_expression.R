#!/usr/bin/env Rscript
# Export aggregated gene x cell_type expression for viewer-side exploration.
# Aggregates only (no cell-level values): mean normalized expression,
# detection rate from raw counts, n_cells, and mean raw counts.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(optparse)
  library(arrow)
  library(jsonlite)
})

option_list <- list(
  make_option("--input",            type = "character", help = "Seurat .rds"),
  make_option("--output_dir",       type = "character", default = "output/export/expression"),
  make_option("--cell_type_column", type = "character", default = "C185_named"),
  make_option("--cell_types",       type = "character", default = ""),
  make_option("--genes",            type = "character", default = ""),
  make_option("--by_dataset",       type = "character", default = "FALSE"),
  make_option("--dataset_column",   type = "character", default = "Dataset"),
  make_option("--datasets",         type = "character", default = ""),
  make_option("--min_cells",        type = "integer",   default = 10L),
  make_option("--assay",            type = "character", default = "RNA"),
  make_option("--value_layer",      type = "character", default = "data"),
  make_option("--counts_layer",     type = "character", default = "counts"),
  make_option("--exclude_col",      type = "character", default = "Author_Exclude"),
  make_option("--exclude_vals",     type = "character", default = "yes"),
  make_option("--drop_zero_rows",   type = "character", default = "TRUE"),
  make_option("--gene_id_type",     type = "character", default = "symbol"),
  make_option("--format",           type = "character", default = "parquet"),
  make_option("--species",          type = "character", default = "mouse")
)
opt <- parse_args(OptionParser(option_list = option_list))

as_bool <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES", "Y")
by_dataset   <- as_bool(opt$by_dataset)
drop_zero    <- as_bool(opt$drop_zero_rows)
stopifnot(nzchar(opt$input), file.exists(opt$input))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

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
  if (is.null(m)) stop("[export] cannot read layer '", layer, "' from assay '", assay, "'")
  m
}

write_table <- function(df, path) {
  if (identical(opt$format, "csv")) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(df, path)
    } else {
      utils::write.csv(df, path, row.names = FALSE)
    }
  } else {
    arrow::write_parquet(df, path, compression = "zstd")
  }
}

message("[export] loading ", opt$input)
obj <- readRDS(opt$input)
if (!inherits(obj, "Seurat")) stop("[export] input is not a Seurat object")
message("[export] ", nrow(obj), " genes x ", ncol(obj), " cells")

md <- obj@meta.data
ct_col <- opt$cell_type_column
if (!ct_col %in% colnames(md)) {
  stop("[export] cell_type_column '", ct_col, "' not found. Available: ",
       paste(colnames(md), collapse = ", "))
}
if (by_dataset && !opt$dataset_column %in% colnames(md)) {
  stop("[export] dataset_column '", opt$dataset_column, "' not found")
}

keep <- !is.na(md[[ct_col]])
excl_vals <- split_arg(opt$exclude_vals)
if (nzchar(opt$exclude_col) && opt$exclude_col %in% colnames(md) && length(excl_vals)) {
  keep <- keep & (is.na(md[[opt$exclude_col]]) | !(md[[opt$exclude_col]] %in% excl_vals))
}
if (by_dataset) keep <- keep & !is.na(md[[opt$dataset_column]])

ct_all <- as.character(md[[ct_col]])
ct_req <- split_arg(opt$cell_types)
if (length(ct_req)) {
  missing_ct <- setdiff(ct_req, unique(ct_all[keep]))
  if (length(missing_ct)) message("[export] WARNING requested cell types absent: ",
                                  paste(missing_ct, collapse = ", "))
  keep <- keep & ct_all %in% ct_req
}

ds_all <- if (opt$dataset_column %in% colnames(md))
  as.character(md[[opt$dataset_column]]) else rep(NA_character_, nrow(md))
ds_req <- split_arg(opt$datasets)
if (by_dataset && length(ds_req)) keep <- keep & ds_all %in% ds_req

genes_req <- split_arg(opt$genes)
genes_use <- if (length(genes_req)) {
  missing_g <- setdiff(genes_req, rownames(obj))
  if (length(missing_g)) message("[export] WARNING requested genes absent: ",
                                 paste(missing_g, collapse = ", "))
  intersect(genes_req, rownames(obj))
} else {
  rownames(obj)
}
if (!length(genes_use)) stop("[export] no genes to export")

kept_idx <- which(keep)
if (!length(kept_idx)) stop("[export] no cells remain after filtering")
ct_kept <- ct_all[kept_idx]
ds_kept <- ds_all[kept_idx]
keys <- if (by_dataset) paste(ct_kept, ds_kept, sep = " || ") else ct_kept

ct_counts <- sort(table(keys), decreasing = TRUE)
keys_u <- names(ct_counts)
n_keys <- length(keys_u)

data_all   <- get_layer(obj, opt$assay, opt$value_layer)
counts_all <- get_layer(obj, opt$assay, opt$counts_layer)

if (inherits(counts_all, "dgCMatrix") && length(counts_all@x)) {
  cx <- counts_all@x
  if (min(cx) < 0) stop("[export] counts contain negative values")
  cs <- if (length(cx) > 1e6) cx[sample.int(length(cx), 1e6)] else cx
  if (any(abs(cs - round(cs)) > 1e-8)) {
    warning("[export] counts are not integer-valued; detection uses counts > 0")
  }
}

genes_n <- length(genes_use)
mean_mat       <- matrix(NA_real_, genes_n, n_keys, dimnames = list(genes_use, keys_u))
mean_counts_mat<- matrix(NA_real_, genes_n, n_keys, dimnames = list(genes_use, keys_u))
pct_mat        <- matrix(NA_real_, genes_n, n_keys, dimnames = list(genes_use, keys_u))
nexpr_mat      <- matrix(NA_real_, genes_n, n_keys, dimnames = list(genes_use, keys_u))
n_cells_vec    <- integer(n_keys); names(n_cells_vec) <- keys_u

message("[export] aggregating ", genes_n, " genes across ", n_keys, " groups")
for (j in seq_along(keys_u)) {
  k   <- keys_u[j]
  idx <- kept_idx[keys == k]
  n   <- length(idx)
  n_cells_vec[j] <- n
  if (n < opt$min_cells) next
  d <- data_all[genes_use, idx, drop = FALSE]
  c <- counts_all[genes_use, idx, drop = FALSE]
  mean_mat[, j]        <- Matrix::rowSums(d) / n
  mean_counts_mat[, j] <- Matrix::rowSums(c) / n
  ne <- Matrix::rowSums(c > 0)
  nexpr_mat[, j] <- ne
  pct_mat[, j]   <- 100 * ne / n
  rm(d, c, ne)
  if (j %% 20 == 0) gc(verbose = FALSE)
}

valid <- n_cells_vec >= opt$min_cells
keys_valid <- keys_u[valid]
if (!length(keys_valid)) stop("[export] no group has >= min_cells (", opt$min_cells, ")")

mean_mat        <- mean_mat[, valid, drop = FALSE]
mean_counts_mat <- mean_counts_mat[, valid, drop = FALSE]
pct_mat         <- pct_mat[, valid, drop = FALSE]
nexpr_mat       <- nexpr_mat[, valid, drop = FALSE]
n_cells_valid   <- n_cells_vec[valid]

ct_valid <- if (by_dataset) sub(" \\|\\| .*$", "", keys_valid) else keys_valid
ds_valid <- if (by_dataset) sub("^.* \\|\\| ", "", keys_valid) else rep(NA_character_, length(keys_valid))

n_cells_input <- length(kept_idx)
datasets_all  <- if (any(!is.na(ds_kept))) sort(unique(ds_kept)) else character(0)
rm(obj, md, data_all, counts_all)
gc(verbose = FALSE)

max_mean <- apply(mean_mat, 1, max, na.rm = TRUE)
max_pct  <- apply(pct_mat, 1, max, na.rm = TRUE)
n_detected <- rowSums(pct_mat > 0, na.rm = TRUE)
top_idx <- max.col(mean_mat, ties.method = "first")
top_ct  <- colnames(mean_mat)[top_idx]
top_mean <- mean_mat[cbind(seq_len(nrow(mean_mat)), top_idx)]

gene_index <- data.frame(
  gene_id                = genes_use,
  gene_name              = genes_use,
  max_mean_expression    = as.numeric(max_mean),
  max_pct_expressing     = as.numeric(max_pct),
  n_cell_types_detected  = as.integer(n_detected),
  top_cell_type          = top_ct,
  top_mean_expression    = as.numeric(top_mean),
  stringsAsFactors = FALSE
)

has_ds <- any(!is.na(ds_kept))
ct_summary <- do.call(rbind, lapply(seq_along(keys_valid), function(j) {
  k <- keys_valid[j]
  n_genes_det <- sum(pct_mat[, j] > 0, na.rm = TRUE)
  if (!has_ds) {
    return(data.frame(
      cell_type = ct_valid[j], n_cells = n_cells_valid[j],
      n_genes_detected = n_genes_det, n_datasets = 0L,
      top_dataset = NA_character_, top_dataset_n = NA_integer_,
      top_dataset_frac = NA_real_, dataset_confined = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  ds_tab <- table(ds_kept[keys == k])
  top_ds <- names(ds_tab)[which.max(ds_tab)]
  data.frame(
    cell_type        = ct_valid[j],
    n_cells          = n_cells_valid[j],
    n_genes_detected = n_genes_det,
    n_datasets       = length(ds_tab),
    top_dataset      = top_ds,
    top_dataset_n    = as.integer(ds_tab[top_ds]),
    top_dataset_frac = as.numeric(ds_tab[top_ds]) / n_cells_valid[j],
    dataset_confined = length(ds_tab) == 1,
    stringsAsFactors = FALSE
  )
}))
if (by_dataset) ct_summary$dataset <- ds_valid

valid_mask <- keys %in% keys_valid
ct_ds <- as.data.frame(
  table(cell_type = ct_kept[valid_mask], dataset = ds_kept[valid_mask]),
  stringsAsFactors = FALSE
)
names(ct_ds)[3] <- "n_cells"
ct_ds <- ct_ds[ct_ds$n_cells > 0, , drop = FALSE]

skipped <- data.frame(
  group    = keys_u[!valid],
  n_cells  = as.integer(n_cells_vec[!valid]),
  reason   = rep("below_min_cells", sum(!valid)),
  stringsAsFactors = FALSE
)
if (by_dataset && nrow(skipped)) {
  skipped$cell_type <- sub(" \\|\\| .*$", "", skipped$group)
  skipped$dataset   <- sub("^.* \\|\\| ", "", skipped$group)
}

long_list <- vector("list", length(keys_valid))
for (j in seq_along(keys_valid)) {
  ne <- nexpr_mat[, j]
  sel <- if (drop_zero) ne > 0 else rep(TRUE, genes_n)
  if (!any(sel)) next
  long_list[[j]] <- data.frame(
    gene_id         = genes_use[sel],
    gene_name       = genes_use[sel],
    cell_type       = ct_valid[j],
    dataset         = ds_valid[j],
    n_cells         = as.integer(n_cells_valid[j]),
    n_expressing    = as.integer(ne[sel]),
    pct_expressing  = as.numeric(pct_mat[sel, j]),
    mean_expression = as.numeric(mean_mat[sel, j]),
    mean_counts     = as.numeric(mean_counts_mat[sel, j]),
    stringsAsFactors = FALSE
  )
}
long <- do.call(rbind, long_list)
if (is.null(long)) {
  long <- data.frame(
    gene_id = character(0), gene_name = character(0), cell_type = character(0),
    dataset = character(0), n_cells = integer(0), n_expressing = integer(0),
    pct_expressing = numeric(0), mean_expression = numeric(0),
    mean_counts = numeric(0), stringsAsFactors = FALSE
  )
}
rm(long_list); gc(verbose = FALSE)

wide_mean <- data.frame(gene_id = genes_use, mean_mat, check.names = FALSE)
wide_pct  <- data.frame(gene_id = genes_use, pct_mat, check.names = FALSE)

files <- list(
  expression_long          = paste0("expression_long.", opt$format),
  expression_wide_mean     = paste0("expression_wide_mean.", opt$format),
  expression_wide_pct      = paste0("expression_wide_pct.", opt$format),
  gene_index               = paste0("gene_index.", opt$format),
  cell_type_summary        = paste0("cell_type_summary.", opt$format),
  cell_type_dataset_counts = paste0("cell_type_dataset_counts.", opt$format),
  skipped_groups           = paste0("skipped_groups.", opt$format),
  manifest                 = "manifest.json",
  columns                  = "COLUMNS.md"
)

write_table(long,       file.path(opt$output_dir, files$expression_long))
write_table(wide_mean,  file.path(opt$output_dir, files$expression_wide_mean))
write_table(wide_pct,   file.path(opt$output_dir, files$expression_wide_pct))
write_table(gene_index, file.path(opt$output_dir, files$gene_index))
write_table(ct_summary, file.path(opt$output_dir, files$cell_type_summary))
write_table(ct_ds,      file.path(opt$output_dir, files$cell_type_dataset_counts))
write_table(skipped,    file.path(opt$output_dir, files$skipped_groups))

manifest <- list(
  schema_version   = "1.0",
  table            = "expression",
  source           = normalizePath(opt$input),
  species          = opt$species,
  gene_id_type     = opt$gene_id_type,
  cell_type_column = ct_col,
  by_dataset       = by_dataset,
  dataset_column   = if (by_dataset) opt$dataset_column else NULL,
  datasets         = datasets_all,
  n_genes          = genes_n,
  n_genes_detected = sum(rowSums(nexpr_mat) > 0, na.rm = TRUE),
  n_cell_types     = length(keys_valid),
  n_cells_input    = n_cells_input,
  n_cells_used     = sum(n_cells_valid),
  min_cells        = opt$min_cells,
  value_layer      = opt$value_layer,
  counts_layer     = opt$counts_layer,
  normalization    = "value_layer as stored in the object (no re-normalization)",
  detection        = "raw counts > 0",
  exclusion        = list(column = opt$exclude_col, values = I(excl_vals)),
  drop_zero_rows   = drop_zero,
  rows             = list(
    expression_long = nrow(long),
    gene_index      = nrow(gene_index),
    cell_type_summary = nrow(ct_summary)
  ),
  files = files
)
jsonlite::write_json(manifest, file.path(opt$output_dir, files$manifest),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")

columns_md <- c(
  "# Expression export columns",
  "",
  paste0("Source: `", normalizePath(opt$input), "`"),
  paste0("Cell type column: `", ct_col, "`"),
  paste0("Gene ID type: `", opt$gene_id_type, "`"),
  "",
  "## expression_long",
  "",
  "| column | description |",
  "|---|---|",
  "| gene_id | gene identifier (see gene_id_type) |",
  "| gene_name | gene symbol |",
  "| cell_type | cell type label |",
  "| dataset | dataset label (NA when by_dataset = FALSE) |",
  "| n_cells | number of cells in this cell type / dataset |",
  "| n_expressing | cells with raw counts > 0 |",
  "| pct_expressing | 100 * n_expressing / n_cells |",
  "| mean_expression | mean of the normalized value layer |",
  "| mean_counts | mean of raw counts |",
  "",
  "## expression_wide_mean / expression_wide_pct",
  "",
  "First column `gene_id`, then one column per cell type. Column names are",
  "the original cell type labels (may contain `:` and spaces).",
  "",
  "## gene_index",
  "",
  "| column | description |",
  "|---|---|",
  "| gene_id, gene_name | gene identifiers |",
  "| max_mean_expression | highest mean expression across cell types |",
  "| max_pct_expressing | highest detection rate across cell types |",
  "| n_cell_types_detected | cell types with pct_expressing > 0 |",
  "| top_cell_type | cell type with the highest mean expression |",
  "| top_mean_expression | mean expression in top_cell_type |",
  "",
  "## cell_type_summary",
  "",
  "| column | description |",
  "|---|---|",
  "| cell_type | cell type label |",
  "| n_cells | cells in this cell type (after filtering) |",
  "| n_genes_detected | genes with pct_expressing > 0 |",
  "| n_datasets | number of datasets contributing |",
  "| top_dataset | dataset with the most cells |",
  "| top_dataset_n | cells from top_dataset |",
  "| top_dataset_frac | top_dataset_n / n_cells |",
  "| dataset_confined | TRUE if only one dataset contributes |",
  "",
  "## Notes",
  "",
  "- Values are aggregated; no cell-level expression is exported.",
  "- `mean_expression` is descriptive, not differential expression.",
  "- Detection is based on raw counts > 0, not normalized values.",
  "- Groups with fewer than `min_cells` cells are listed in `skipped_groups`.",
  "- With `drop_zero_rows = TRUE`, long-format rows with n_expressing = 0 are omitted;",
  "  `gene_index` still contains every gene in the input (use it for gene search).",
  ""
)
writeLines(columns_md, file.path(opt$output_dir, files$columns))

invisible(file.create(file.path(opt$output_dir, "expression_complete.flag")))
message("[export] done: ", opt$output_dir, " | long rows: ", nrow(long),
        " | cell types: ", length(keys_valid))
