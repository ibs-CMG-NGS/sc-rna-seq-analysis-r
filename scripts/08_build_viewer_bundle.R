#!/usr/bin/env Rscript
# Assemble a cmg-seqviewer handoff bundle from an expression export.
# Layout follows the viewer team's 2026-10-20 reply:
#   <bundle>/seqviewer_manifest.json
#   <bundle>/datasets/<prefix>_wide_mean.parquet
#   <bundle>/datasets/<prefix>_wide_pct.parquet
#   <bundle>/datasets/hypomap_gene_index.parquet
#   <bundle>/datasets/hypomap_cell_type_summary.parquet
# `file_path` in the manifest is a bare filename (no directories).

suppressPackageStartupMessages({
  library(arrow)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--export_dir",   type = "character", help = "07_export_expression.R output dir"),
  make_option("--bundle_dir",   type = "character", help = "viewer bundle root"),
  make_option("--prefix",       type = "character", default = "hypomap"),
  make_option("--alias",        type = "character", default = ""),
  make_option("--dataset_type", type = "character", default = "expression_summary"),
  make_option("--primary",      type = "character", default = "wide"),
  make_option("--cell_type_level", type = "character", default = ""),
  make_option("--organism",     type = "character", default = "Mus musculus"),
  make_option("--source",       type = "character", default = "HypoMap"),
  make_option("--license",      type = "character", default = "CC BY 4.0"),
  make_option("--citation",     type = "character", default = ""),
  make_option("--pct_scale",    type = "character", default = "0-100"),
  make_option("--gene_id_type", type = "character", default = "symbol"),
  make_option("--normalization", type = "character", default = "log1p_normalized"),
  make_option("--detection",    type = "character", default = "count > 0"),
  make_option("--gene_index_name", type = "character", default = "hypomap_gene_index.parquet"),
  make_option("--summary_name", type = "character", default = "hypomap_cell_type_summary.parquet"),
  make_option("--include_long", type = "character", default = "FALSE"),
  make_option("--long_name",    type = "character", default = ""),
  make_option("--fallback_multi_group", type = "character", default = "FALSE",
              help = "Deprecated: viewer >=1.4 supports expression_summary [default: %default]"),
  make_option("--tags",         type = "character", default = ""),
  make_option("--notes",        type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

as_bool <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES", "Y")
stopifnot(nzchar(opt$export_dir), dir.exists(opt$export_dir),
          nzchar(opt$bundle_dir), nzchar(opt$prefix),
          opt$primary %in% c("wide", "long"))
include_long <- as_bool(opt$include_long)
fallback_mg  <- as_bool(opt$fallback_multi_group)
if (!nzchar(opt$alias)) opt$alias <- opt$prefix
if (!nzchar(opt$long_name)) opt$long_name <- paste0(opt$prefix, "_long.parquet")

ds_dir <- file.path(opt$bundle_dir, "datasets")
dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

read_export <- function(name, required = TRUE) {
  p <- file.path(opt$export_dir, name)
  if (!file.exists(p)) {
    if (required) stop("[bundle] missing ", p)
    return(NULL)
  }
  arrow::read_parquet(p)
}

manifest_in <- fromJSON(file.path(opt$export_dir, "manifest.json"))
wide_mean <- read_export("expression_wide_mean.parquet")
wide_pct  <- read_export("expression_wide_pct.parquet")
gene_idx  <- read_export("gene_index.parquet")
ct_sum    <- read_export("cell_type_summary.parquet")

if (!"gene_name" %in% names(wide_mean)) wide_mean$gene_name <- wide_mean$gene_id
if (!"gene_name" %in% names(wide_pct))  wide_pct$gene_name  <- wide_pct$gene_id
reorder_cols <- function(df) {
  df[, c("gene_id", "gene_name",
         setdiff(names(df), c("gene_id", "gene_name"))), drop = FALSE]
}
wide_mean <- reorder_cols(wide_mean)
wide_pct  <- reorder_cols(wide_pct)
if ("max_mean_expression" %in% names(gene_idx))
  names(gene_idx)[names(gene_idx) == "max_mean_expression"] <- "max_mean"
if ("max_pct_expressing" %in% names(gene_idx))
  names(gene_idx)[names(gene_idx) == "max_pct_expressing"] <- "max_pct"
if (!"gene_name" %in% names(gene_idx)) gene_idx$gene_name <- gene_idx$gene_id

rownames(wide_mean) <- NULL
rownames(wide_pct)  <- NULL

write_table <- function(df, name) {
  path <- file.path(ds_dir, name)
  arrow::write_parquet(df, path, compression = "zstd")
  path
}
mean_path <- write_table(wide_mean, paste0(opt$prefix, "_wide_mean.parquet"))
pct_path  <- write_table(wide_pct,  paste0(opt$prefix, "_wide_pct.parquet"))
gi_path   <- write_table(gene_idx,  opt$gene_index_name)
sum_path  <- write_table(ct_sum,    opt$summary_name)

uuid_like <- function(x) {
  f <- tempfile(); writeLines(x, f)
  h <- substr(unname(tools::md5sum(f)), 1, 32); unlink(f)
  paste(substr(h, 1, 8), substr(h, 9, 12), substr(h, 13, 16),
        substr(h, 17, 20), substr(h, 21, 32), sep = "-")
}

tags <- trimws(strsplit(opt$tags, ",")[[1]])
tags <- tags[nzchar(tags)]
if (!length(tags)) tags <- c(opt$source, opt$cell_type_level, "expression")

base_notes <- if (nzchar(opt$notes)) opt$notes else
  paste0(opt$source, " ", opt$cell_type_level, "; mean normalized expression")
notes <- paste0(base_notes, " | n_cells_used=", manifest_in$n_cells_used,
                ", min_cells=", manifest_in$min_cells,
                ", export_schema=", manifest_in$schema_version)

make_schema <- function(layout) {
  companions <- if (layout == "wide") {
    list(pct = basename(pct_path), gene_index = basename(gi_path),
         cell_type_summary = basename(sum_path))
  } else {
    list(gene_index = basename(gi_path), cell_type_summary = basename(sum_path))
  }
  list(
    schema_version        = "1.0",
    layout                = layout,
    value_type            = if (layout == "wide") "mean_expression" else "long",
    cell_type_level       = opt$cell_type_level,
    gene_id_type          = opt$gene_id_type,
    normalization         = opt$normalization,
    detection             = opt$detection,
    pct_scale             = opt$pct_scale,
    companion_files       = companions,
    value_columns         = list(
      mean_expression = "mean_expression",
      pct_expressing  = "pct_expressing",
      n_cells         = "n_cells",
      n_expressing    = "n_expressing",
      mean_counts     = "mean_counts"
    ),
    source                = opt$source,
    license               = opt$license,
    citation              = opt$citation,
    export_schema_version = manifest_in$schema_version,
    cell_type_column      = manifest_in$cell_type_column,
    n_cells_used          = manifest_in$n_cells_used,
    min_cells             = manifest_in$min_cells
  )
}

make_item <- function(df, path, layout, alias_suffix = "", item_notes = notes) {
  n_genes <- if (layout == "wide") nrow(df) else length(unique(df$gene_id))
  n_cts   <- if (layout == "wide") ncol(df) - 2L else length(unique(df$cell_type))
  list(
    dataset_id        = uuid_like(basename(path)),
    alias             = paste0(opt$alias, alias_suffix),
    file_path         = basename(path),
    original_filename = basename(path),
    dataset_type      = opt$dataset_type,
    organism          = opt$organism,
    row_count         = nrow(df),
    gene_count        = n_genes,
    cell_type_count   = n_cts,
    notes             = item_notes,
    tags              = as.list(tags),
    expression_schema = make_schema(layout)
  )
}

long_path <- NULL
if (opt$primary == "long" || include_long) {
  long <- read_export("expression_long.parquet")
  long_path <- write_table(long, opt$long_name)
}

if (opt$primary == "long") {
  main_item <- make_item(long, long_path, "long")
} else {
  main_item <- make_item(wide_mean, mean_path, "wide")
}
manifest <- list(schema_version = "1.0", datasets = list(main_item))

if (include_long && opt$primary == "wide") {
  long_item <- make_item(
    long, long_path, "long", " (long)",
    paste0(notes, " | long layout (per cell_type x dataset)")
  )
  manifest$datasets <- c(manifest$datasets, list(long_item))
}

jsonlite::write_json(manifest, file.path(opt$bundle_dir, "seqviewer_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")

if (fallback_mg) {
  mg <- wide_mean[, c("gene_id", setdiff(names(wide_mean), c("gene_id", "gene_name"))),
                  drop = FALSE]
  mg_file <- paste0(opt$prefix, "_wide_mean_multi_group.parquet")
  mg_path <- write_table(mg, mg_file)
  mg_item <- list(
    dataset_id        = uuid_like(basename(mg_path)),
    alias             = paste0(opt$alias, " (multi_group test)"),
    file_path         = basename(mg_path),
    original_filename = basename(mg_path),
    dataset_type      = "multi_group",
    organism          = opt$organism,
    row_count         = nrow(mg),
    gene_count        = nrow(mg),
    cell_type_count   = ncol(mg) - 1L,
    notes             = paste0(notes, " | deprecated multi_group fallback (viewer >=1.4 supports expression_summary)"),
    tags              = as.list(c(tags, "multi_group-test")),
    expression_schema = make_schema("wide")
  )
  jsonlite::write_json(list(schema_version = "1.0", datasets = list(mg_item)),
                       file.path(opt$bundle_dir, "seqviewer_manifest.multi_group.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null")
}

readme <- c(
  paste0("# ", opt$alias),
  "",
  "Handoff bundle for cmg-seqviewer (expression summary).",
  "",
  "## Files",
  "",
  paste0("- `datasets/", basename(mean_path), "`: wide mean expression"),
  paste0("- `datasets/", basename(pct_path), "`: wide pct_expressing (", opt$pct_scale, ")"),
  paste0("- `datasets/", basename(gi_path), "`: gene index"),
  paste0("- `datasets/", basename(sum_path), "`: cell type summary"),
  if (!is.null(long_path)) paste0("- `datasets/", basename(long_path), "`: long layout (Phase 2)"),
  "",
  "## Manifest",
  "",
  paste0("- `seqviewer_manifest.json`: final entry (`dataset_type: ", opt$dataset_type, "`)."),
  if (fallback_mg) c(
    "- `seqviewer_manifest.multi_group.json`: deprecated Phase 0 fallback entry.",
    "  viewer >=1.4 supports `expression_summary`; use the main manifest instead."
  ),
  "",
  paste0("pct scale: ", opt$pct_scale, " | gene id: ", opt$gene_id_type,
         " | n_cells_used: ", manifest_in$n_cells_used),
  ""
)
writeLines(readme, file.path(opt$bundle_dir, "README.md"))

message("[bundle] ", opt$bundle_dir, " | primary: ", opt$primary,
        " | files: ", length(list.files(ds_dir)))
