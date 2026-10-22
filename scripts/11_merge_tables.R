#!/usr/bin/env Rscript
# 11_merge_tables.R
# Merge canonical analysis and expression CSV manifests into one table root:
#   <tables_dir>/manifest.json
#   <tables_dir>/index.tsv
#   <tables_dir>/README.md

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

opt_list <- list(
  make_option("--tables_dir", type = "character", help = "Canonical table root"),
  make_option("--analysis_manifest", type = "character", default = ""),
  make_option("--expression_manifest", type = "character", default = ""),
  make_option("--table_schema_version", type = "character", default = "1.0")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (!nzchar(opt$tables_dir)) stop("--tables_dir is required", call. = FALSE)
if (!nzchar(opt$analysis_manifest)) {
  opt$analysis_manifest <- file.path(opt$tables_dir, "analysis", "manifest.json")
}
if (!nzchar(opt$expression_manifest)) {
  opt$expression_manifest <- file.path(opt$tables_dir, "expression", "manifest.json")
}

as_scalar <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  v <- unlist(x, use.names = FALSE)
  if (!length(v)) return(default)
  as.character(v[[1]])
}

as_int <- function(x) suppressWarnings(as.integer(as_scalar(x, NA_character_)))

read_source <- function(path, group) {
  if (!file.exists(path)) {
    return(list(entries = list(), source = NULL))
  }
  m <- fromJSON(path, simplifyVector = FALSE)
  items <- m$datasets
  if (is.null(items) || !length(items)) {
    return(list(entries = list(), source = m))
  }

  entries <- lapply(items, function(item) {
    rel <- as_scalar(item$file_path, "")
    if (!nzchar(rel)) stop("[merge_tables] missing file_path in ", path, call. = FALSE)
    rel <- gsub("\\\\", "/", rel)
    full <- file.path(opt$tables_dir, group, rel)
    if (!file.exists(full)) {
      stop("[merge_tables] missing table for ", group, ": ", full, call. = FALSE)
    }
    item$file_path <- paste0(group, "/", rel)
    item$source_group <- group
    item$table_name <- as_scalar(item$table_name, sub("\\.csv$", "", basename(rel)))
    item$row_count <- as_int(item$row_count)
    item$column_count <- as_int(item$column_count)
    item
  })
  list(entries = entries, source = m)
}

analysis <- read_source(opt$analysis_manifest, "analysis")
expression <- read_source(opt$expression_manifest, "expression")
entries <- c(analysis$entries, expression$entries)

if (length(entries)) {
  ord <- order(
    vapply(entries, function(x) as_scalar(x$source_group, ""), character(1)),
    vapply(entries, function(x) as_scalar(x$dataset_type, ""), character(1)),
    vapply(entries, function(x) as_scalar(x$table_name, ""), character(1))
  )
  entries <- entries[ord]
}

sources <- list()
if (!is.null(analysis$source)) {
  sources[[length(sources) + 1]] <- list(
    group = "analysis",
    manifest = "analysis/manifest.json",
    source = as_scalar(analysis$source$source, "06_export_seqviewer.R"),
    table_schema_version = as_scalar(analysis$source$table_schema_version, NA_character_),
    n_datasets = length(analysis$entries)
  )
}
if (!is.null(expression$source)) {
  sources[[length(sources) + 1]] <- list(
    group = "expression",
    manifest = "expression/manifest.json",
    source = as_scalar(expression$source$source, "07_export_expression.R"),
    table_schema_version = as_scalar(expression$source$table_schema_version, NA_character_),
    n_datasets = length(expression$entries)
  )
}

dir.create(opt$tables_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- list(
  table_schema_version = opt$table_schema_version,
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = entries,
  sources = sources
)
jsonlite::write_json(manifest, file.path(opt$tables_dir, "manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")

index <- if (length(entries)) {
  do.call(rbind, lapply(entries, function(item) {
    data.frame(
      file         = as_scalar(item$file_path, ""),
      source_group = as_scalar(item$source_group, ""),
      dataset_type = as_scalar(item$dataset_type, ""),
      table_name   = as_scalar(item$table_name, ""),
      layout       = as_scalar(item$layout, ""),
      cell_type    = as_scalar(item$cell_type, ""),
      condition    = as_scalar(item$experiment_condition, ""),
      direction    = as_scalar(item$direction, ""),
      ontology     = as_scalar(item$ontology, ""),
      rows         = as_int(item$row_count),
      columns      = as_int(item$column_count),
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    file = character(0), source_group = character(0), dataset_type = character(0),
    table_name = character(0), layout = character(0), cell_type = character(0),
    condition = character(0), direction = character(0), ontology = character(0),
    rows = integer(0), columns = integer(0), stringsAsFactors = FALSE
  )
}
rownames(index) <- NULL
write.table(index, file.path(opt$tables_dir, "index.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

counts <- table(index$dataset_type)
count_lines <- if (length(counts)) {
  paste0("- `", names(counts), "`: ", as.integer(counts), " table(s)")
} else {
  "- (no tables)"
}
readme <- c(
  "# Canonical tables",
  "",
  paste0("Schema: `", opt$table_schema_version, "`"),
  "",
  "This directory is the source-of-truth table layer. Each CSV is a flat table",
  "suitable for Excel, R, Python, or the cmg-seqviewer. Viewer-specific parquet",
  "files are generated separately under `export/seqviewer`.",
  "",
  "## Layout",
  "",
  "- `manifest.json`: merged table manifest (paths are relative to this directory).",
  "- `index.tsv`: one row per table for quick inspection.",
  "- `analysis/datasets/*.csv`: DE, GO/KEGG, composition, and pseudobulk counts.",
  "- `expression/datasets/*.csv`: gene x cell_type expression summaries.",
  "",
  "Excel workbooks are written to `export/excel/` when `table_export.excel` is true.",
  "Long tables above Excel's row limit (1,048,575 data rows) remain CSV-only and",
  "are listed as skipped in `export/excel/excel_manifest.json`.",
  "",
  "## Tables",
  "",
  count_lines,
  ""
)
writeLines(readme, file.path(opt$tables_dir, "README.md"))

message("[merge_tables] ", length(analysis$entries), " analysis + ",
        length(expression$entries), " expression table(s) -> ", opt$tables_dir)
