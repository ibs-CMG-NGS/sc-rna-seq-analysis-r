#!/usr/bin/env Rscript
# 10_merge_viewer_bundle.R
# Merge analysis (DE/GO/composition/MG) and expression-summary staging bundles
# into one cmg-seqviewer Import Folder.
#
# Inputs are produced by:
#   - scripts/06_export_seqviewer.R  (analysis staging)
#   - scripts/08_build_viewer_bundle.R (expression staging)
#   - scripts/13_export_distribution.R (pseudobulk distribution staging)
#
# Output layout:
#   <output_dir>/seqviewer_manifest.json
#   <output_dir>/datasets/*.parquet
#   <output_dir>/export_index.tsv
#   <output_dir>/COLUMNS.md
#   <output_dir>/export_complete.flag

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

opt_list <- list(
  make_option("--analysis_dir",     type = "character", default = ""),
  make_option("--expression_dir",   type = "character", default = ""),
  make_option("--distribution_dir", type = "character", default = ""),
  make_option("--output_dir",       type = "character", help = "Final Import Folder"),
  make_option("--schema_version",   type = "character", default = "1.0")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (!nzchar(opt$output_dir)) stop("--output_dir is required", call. = FALSE)

datasets_dir <- file.path(opt$output_dir, "datasets")
dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(datasets_dir, pattern = "\\.parquet$", full.names = TRUE))
unlink(file.path(opt$output_dir, "seqviewer_manifest.multi_group.json"))

as_scalar <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  v <- unlist(x, use.names = FALSE)
  if (!length(v)) return(default)
  as.character(v[[1]])
}

read_manifest <- function(dir) {
  if (!nzchar(dir) || !dir.exists(dir)) return(list())
  p <- file.path(dir, "seqviewer_manifest.json")
  if (!file.exists(p)) return(list())
  m <- fromJSON(p, simplifyVector = FALSE)
  if (is.null(m$datasets) || !length(m$datasets)) return(list())
  m$datasets
}

copy_one <- function(fp, src_dir, label, seen_files, duplicate_is_error) {
  fp <- basename(as_scalar(fp, ""))
  if (!nzchar(fp)) stop("[merge] missing file_path in ", label, " manifest", call. = FALSE)
  if (exists(fp, envir = seen_files, inherits = FALSE)) {
    if (duplicate_is_error) {
      stop("[merge] duplicate file_path across bundles: ", fp, call. = FALSE)
    }
    return(invisible(NULL))
  }
  src <- file.path(src_dir, "datasets", fp)
  if (!file.exists(src)) src <- file.path(src_dir, fp)
  if (!file.exists(src)) {
    stop("[merge] missing parquet for ", label, ": ", fp, call. = FALSE)
  }
  dst <- file.path(datasets_dir, fp)
  if (!file.copy(src, dst, overwrite = FALSE)) {
    stop("[merge] failed to copy ", src, " -> ", dst, call. = FALSE)
  }
  assign(fp, label, envir = seen_files)
  invisible(NULL)
}

copy_items <- function(items, src_dir, label, seen_files, seen_ids) {
  if (!length(items)) return(invisible(NULL))
  for (item in items) {
    fp <- as_scalar(item$file_path, "")
    copy_one(fp, src_dir, label, seen_files, duplicate_is_error = TRUE)

    companions <- item$expression_schema$companion_files
    if (is.null(companions)) companions <- item$distribution_schema$companion_files
    if (!is.null(companions) && length(companions)) {
      for (cfp in unlist(companions, use.names = FALSE)) {
        copy_one(cfp, src_dir, label, seen_files, duplicate_is_error = FALSE)
      }
    }

    id <- as_scalar(item$dataset_id, "")
    if (nzchar(id)) {
      if (exists(id, envir = seen_ids, inherits = FALSE)) {
        stop("[merge] duplicate dataset_id across bundles: ", id, call. = FALSE)
      }
      assign(id, label, envir = seen_ids)
    }
  }
}

analysis_items     <- read_manifest(opt$analysis_dir)
expression_items   <- read_manifest(opt$expression_dir)
distribution_items <- read_manifest(opt$distribution_dir)
entries <- c(analysis_items, expression_items, distribution_items)
if (!length(entries)) stop("[merge] no datasets to merge", call. = FALSE)

seen_files <- new.env(parent = emptyenv())
seen_ids   <- new.env(parent = emptyenv())
copy_items(analysis_items,     opt$analysis_dir,     "analysis",     seen_files, seen_ids)
copy_items(expression_items,   opt$expression_dir,   "expression",   seen_files, seen_ids)
copy_items(distribution_items, opt$distribution_dir, "distribution", seen_files, seen_ids)

write_json(
  list(schema_version = opt$schema_version, datasets = entries),
  file.path(opt$output_dir, "seqviewer_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, null = "null"
)

rows <- lapply(entries, function(item) {
  data.frame(
    file         = as_scalar(item$file_path, ""),
    dataset_type = as_scalar(item$dataset_type, ""),
    cell_type    = as_scalar(item$cell_type, ""),
    condition    = as_scalar(item$experiment_condition, ""),
    direction    = as_scalar(item$direction, ""),
    ontology     = as_scalar(item$ontology, ""),
    rows         = suppressWarnings(as.integer(as_scalar(item$row_count, NA_character_))),
    stringsAsFactors = FALSE
  )
})
index <- do.call(rbind, rows)
index <- index[order(index$dataset_type, index$cell_type, index$condition, index$file,
                     na.last = TRUE), , drop = FALSE]
rownames(index) <- NULL
write.table(index, file.path(opt$output_dir, "export_index.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

columns_path <- file.path(opt$analysis_dir, "COLUMNS.md")
columns_md <- if (nzchar(opt$analysis_dir) && file.exists(columns_path)) {
  readLines(columns_path, warn = FALSE)
} else {
  c(
    "# cmg-seqviewer handoff: column contract",
    "",
    paste0("schema_version: ", opt$schema_version),
    "",
    "Producer: scRNA-seq 3rd analysis pipeline. Consumer: cmg-seqviewer.",
    "The bundle root is an Import Folder: `seqviewer_manifest.json` + `datasets/*.parquet`."
  )
}

expr_items <- Filter(function(x) identical(as_scalar(x$dataset_type), "expression_summary"), entries)
if (length(expr_items)) {
  expr_layouts <- unique(vapply(expr_items, function(x) {
    layout <- as_scalar(x$expression_schema$layout)
    if (is.na(layout) || !nzchar(layout)) "wide" else layout
  }, character(1)))
  if ("wide" %in% expr_layouts) {
    columns_md <- c(
      columns_md,
      "",
      "## expression_summary (wide gene x cell_type)",
      "",
      "| column | notes |",
      "| --- | --- |",
      "| gene_id, gene_name | gene identifier and symbol |",
      "| <cell_type> | mean normalized expression (numeric) |",
      "",
      "The wide mean table is the primary dataset when `layout: wide`. Companion files",
      "(percent expressing, gene index, cell-type summary) are listed in",
      "`expression_schema.companion_files` for each expression dataset in",
      "`seqviewer_manifest.json`."
    )
  }
  if ("long" %in% expr_layouts) {
    columns_md <- c(
      columns_md,
      "",
      "## expression_summary (long gene x cell_type x dataset)",
      "",
      "| column | notes |",
      "| --- | --- |",
      "| gene_id, gene_name | gene identifier and symbol |",
      "| cell_type | cell type label |",
      "| dataset | source dataset (when `by_dataset: true`) |",
      "| n_cells, n_expressing | cells in group, cells with nonzero expression |",
      "| pct_expressing | percent expressing, scale from `expression_schema.pct_scale` |",
      "| mean_expression | mean normalized expression |",
      "| mean_counts | mean raw counts |",
      "",
      "The long table is the primary dataset when `layout: long`. Companion files",
      "(gene index, cell-type summary) are listed in",
      "`expression_schema.companion_files` for each expression dataset in",
      "`seqviewer_manifest.json`."
    )
  }
}

dist_items <- Filter(function(x) identical(as_scalar(x$dataset_type), "expression_distribution"), entries)
if (length(dist_items)) {
  columns_md <- c(
    columns_md,
    "",
    "## expression_distribution (per-sample pseudobulk)",
    "",
    "One row = gene x cell_type x sample. Companion files are listed in",
    "`distribution_schema.companion_files`.",
    "",
    "| column | notes |",
    "| --- | --- |",
    "| gene_id, gene_name | gene symbol |",
    "| cell_type | cell type label |",
    "| sample | biological sample/replicate ID |",
    "| dataset, condition | study and experimental condition |",
    "| n_cells, n_expressing | cells in group, cells with raw counts > 0 |",
    "| pct_expressing | percent expressing (0-100) |",
    "| mean_expression | mean normalized value layer |",
    "| mean_counts | mean raw counts |"
  )
}
writeLines(columns_md, file.path(opt$output_dir, "COLUMNS.md"))

invisible(file.create(file.path(opt$output_dir, "export_complete.flag")))

message("[merge] ", length(analysis_items), " analysis + ",
        length(expression_items), " expression + ",
        length(distribution_items), " distribution dataset(s) -> ", opt$output_dir)
