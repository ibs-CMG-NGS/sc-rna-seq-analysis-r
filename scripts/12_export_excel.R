#!/usr/bin/env Rscript
# 12_export_excel.R
# Build one Excel workbook per dataset_type from the canonical CSV table layer.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

opt_list <- list(
  make_option("--tables_dir", type = "character", help = "Canonical table root"),
  make_option("--output_dir", type = "character", default = "",
              help = "Excel output directory [default: <tables_dir>/../export/excel]"),
  make_option("--manifest", type = "character", default = ""),
  make_option("--max_rows", type = "integer", default = 1048575L,
              help = "Maximum data rows per sheet (Excel limit minus header)"),
  make_option("--max_cols", type = "integer", default = 16383L,
              help = "Maximum data columns per sheet (Excel limit)")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (!nzchar(opt$tables_dir)) stop("--tables_dir is required", call. = FALSE)
if (!nzchar(opt$manifest)) opt$manifest <- file.path(opt$tables_dir, "manifest.json")
if (!nzchar(opt$output_dir)) {
  opt$output_dir <- file.path(dirname(normalizePath(opt$tables_dir, mustWork = FALSE)),
                              "export", "excel")
}
if (!file.exists(opt$manifest)) stop("missing manifest: ", opt$manifest, call. = FALSE)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("R package 'openxlsx' is required", call. = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("R package 'data.table' is required", call. = FALSE)
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fromJSON(opt$manifest, simplifyVector = FALSE)
entries <- manifest$datasets
if (is.null(entries)) entries <- list()

as_scalar <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  v <- unlist(x, use.names = FALSE)
  if (!length(v)) return(default)
  as.character(v[[1]])
}
as_int <- function(x) suppressWarnings(as.integer(as_scalar(x, NA_character_)))

safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
make_unique <- function(x, max_chars = 31L) {
  x <- substr(safe(x), 1L, max_chars)
  x[!nzchar(x)] <- "table"
  out <- character(length(x))
  used <- character(0)
  for (i in seq_along(x)) {
    candidate <- x[i]
    n <- 1L
    while (tolower(candidate) %in% used) {
      suffix <- paste0("_", n)
      candidate <- paste0(substr(x[i], 1L, max_chars - nchar(suffix)), suffix)
      n <- n + 1L
    }
    used <- c(used, tolower(candidate))
    out[i] <- candidate
  }
  out
}

types <- vapply(entries, function(x) as_scalar(x$dataset_type, "table"), character(1))
workbooks <- list()
skipped <- list()

for (type in unique(types)) {
  idx <- which(types == type)
  items <- entries[idx]
  wb <- openxlsx::createWorkbook()
  sheets <- list()
  item_names <- vapply(items, function(x) {
    nm <- as_scalar(x$table_name, "")
    if (nzchar(nm)) nm else as_scalar(x$dataset_type, "table")
  }, character(1))
  sheet_names <- make_unique(item_names)

  for (k in seq_along(items)) {
    item <- items[[k]]
    table_name <- item_names[k]
    sheet <- sheet_names[k]
    rel <- as_scalar(item$file_path, "")
    if (!nzchar(rel)) {
      skipped[[length(skipped) + 1]] <- list(
        table_name = table_name, dataset_type = type, reason = "missing_file_path"
      )
      next
    }
    csv_path <- file.path(opt$tables_dir, rel)
    if (!file.exists(csv_path)) {
      skipped[[length(skipped) + 1]] <- list(
        table_name = table_name, dataset_type = type, file_path = rel,
        reason = "missing_csv"
      )
      next
    }

    df <- data.table::fread(csv_path, data.table = FALSE)
    rows <- nrow(df)
    cols <- ncol(df)
    if (rows > opt$max_rows) {
      skipped[[length(skipped) + 1]] <- list(
        table_name = table_name, dataset_type = type, file_path = rel,
        rows = rows, columns = cols,
        reason = paste0("exceeds_excel_row_limit_", opt$max_rows)
      )
      next
    }
    if (cols > opt$max_cols) {
      skipped[[length(skipped) + 1]] <- list(
        table_name = table_name, dataset_type = type, file_path = rel,
        rows = rows, columns = cols,
        reason = paste0("exceeds_excel_column_limit_", opt$max_cols)
      )
      next
    }

    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df, keepNA = TRUE, na.string = "")
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    if (rows > 0L && cols > 0L) {
      openxlsx::addFilter(wb, sheet, row = 1L, cols = seq_len(cols))
    }
    sheets[[length(sheets) + 1]] <- list(
      sheet = sheet,
      table_name = table_name,
      file_path = rel,
      rows = rows,
      columns = cols
    )
  }

  if (!length(sheets)) {
    skipped[[length(skipped) + 1]] <- list(
      dataset_type = type, reason = "no_exportable_sheets"
    )
    next
  }

  wb_name <- paste0(make_unique(type)[1], ".xlsx")
  wb_path <- file.path(opt$output_dir, wb_name)
  openxlsx::saveWorkbook(wb, wb_path, overwrite = TRUE)
  workbooks[[length(workbooks) + 1]] <- list(
    file = wb_name,
    dataset_type = type,
    sheets = sheets,
    n_sheets = length(sheets)
  )
}

excel_manifest <- list(
  excel_schema_version = "1.0",
  table_schema_version = as_scalar(manifest$table_schema_version, NA_character_),
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  limits = list(max_rows = opt$max_rows, max_cols = opt$max_cols),
  workbooks = workbooks,
  skipped = skipped
)
jsonlite::write_json(excel_manifest, file.path(opt$output_dir, "excel_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")

wb_lines <- if (length(workbooks)) {
  vapply(workbooks, function(wb) {
    paste0("- `", wb$file, "`: ", wb$n_sheets, " sheet(s) (",
           as_scalar(wb$dataset_type, ""), ")")
  }, character(1))
} else {
  "- (no workbooks)"
}
skip_lines <- if (length(skipped)) {
  vapply(skipped, function(x) {
    rows <- as_int(x$rows)
    cols <- as_int(x$columns)
    dims <- if (!is.na(rows) && !is.na(cols)) paste0(rows, " x ", cols) else ""
    paste0("- `", as_scalar(x$table_name, as_scalar(x$dataset_type, "table")), "`",
           if (nzchar(dims)) paste0(" (", dims, ")") else "",
           ": ", as_scalar(x$reason, "skipped"))
  }, character(1))
} else {
  "- (none)"
}
readme <- c(
  "# Excel export",
  "",
  paste0("Workbook set built from the canonical CSV tables (`",
         as_scalar(manifest$table_schema_version, "1.0"), "` schema)."),
  "",
  "## Workbooks",
  "",
  wb_lines,
  "",
  "## Skipped tables",
  "",
  skip_lines,
  "",
  paste0("Excel limits: ", opt$max_rows, " data rows and ", opt$max_cols,
         " data columns per sheet. Larger tables remain CSV-only; see ",
         "`excel_manifest.json` for details."),
  ""
)
writeLines(readme, file.path(opt$output_dir, "README.md"))

message("[export_excel] ", length(workbooks), " workbook(s), ",
        length(skipped), " skipped -> ", opt$output_dir)
