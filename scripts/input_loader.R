#!/usr/bin/env Rscript
# input_loader.R
# Load and validate Seurat (.rds) or AnnData (.h5ad) input for the pipeline.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
})

# ── CLI args ──────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",         type = "character", help = "Path to input file (.rds or .h5ad)"),
  make_option("--format",        type = "character", default = "seurat",
              help = "Input format: seurat | anndata [default: %default]"),
  make_option("--cell_type_col", type = "character", default = "cell_type",
              help = "Metadata column name for cell type [default: %default]"),
  make_option("--sample_col",    type = "character", default = "sample_id",
              help = "Metadata column name for sample [default: %default]"),
  make_option("--condition_col", type = "character", default = "condition",
              help = "Metadata column name for condition [default: %default]"),
  make_option("--cell_types",    type = "character", default = "null",
              help = "Comma-separated cell types to include; 'null' = all [default: %default]"),
  make_option("--min_cells",     type = "integer",   default = 10L,
              help = "Minimum cells per sample×cell_type group [default: %default]"),
  make_option("--output_rds",    type = "character", help = "Path to save validated Seurat .rds")
)

opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$input),      "--input is required")
stopifnot(!is.null(opt$output_rds), "--output_rds is required")

# ── Load ──────────────────────────────────────────────────────────────────────
message("[input_loader] Loading input: ", opt$input)

if (opt$format == "seurat") {
  obj <- readRDS(opt$input)
  if (!inherits(obj, "Seurat")) {
    stop("Object loaded from .rds is not a Seurat object.")
  }

} else if (opt$format == "anndata") {
  suppressPackageStartupMessages(library(reticulate))
  ad <- reticulate::import("anndata")
  sc <- reticulate::import("scanpy")

  adata <- ad$read_h5ad(opt$input)

  # Convert to Seurat via counts matrix
  counts <- t(as.matrix(adata$X))
  rownames(counts) <- rownames(adata$var)
  colnames(counts) <- rownames(adata$obs)

  obj <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(adata$obs))
  message("[input_loader] Converted AnnData -> Seurat: ", ncol(obj), " cells, ", nrow(obj), " genes")

} else {
  stop("Unknown input format: ", opt$format, ". Use 'seurat' or 'anndata'.")
}

# ── Validate required metadata columns ───────────────────────────────────────
required_cols <- c(opt$cell_type_col, opt$sample_col, opt$condition_col)
missing_cols  <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(colnames(obj@meta.data), collapse = ", "))
}

# ── Filter cell types (if requested) ─────────────────────────────────────────
if (!is.null(opt$cell_types) && tolower(opt$cell_types) != "null" && nchar(trimws(opt$cell_types)) > 0) {
  keep_types <- trimws(strsplit(opt$cell_types, ",")[[1]])
  all_types  <- unique(obj@meta.data[[opt$cell_type_col]])
  bad_types  <- setdiff(keep_types, all_types)
  if (length(bad_types) > 0) {
    warning("Requested cell types not found in data: ", paste(bad_types, collapse = ", "))
  }
  obj <- subset(obj, cells = which(obj@meta.data[[opt$cell_type_col]] %in% keep_types))
  message("[input_loader] Filtered to ", length(keep_types), " cell type(s): ", ncol(obj), " cells remain")
}

# ── Check min cells per group ─────────────────────────────────────────────────
group_counts <- table(
  obj@meta.data[[opt$cell_type_col]],
  obj@meta.data[[opt$sample_col]]
)
small_groups <- which(group_counts < opt$min_cells & group_counts > 0, arr.ind = TRUE)
if (nrow(small_groups) > 0) {
  warning(
    "[input_loader] ", nrow(small_groups), " cell_type × sample group(s) have < ",
    opt$min_cells, " cells (will be excluded from pseudobulk DE):\n",
    paste(apply(small_groups, 1, function(r) {
      paste0("  ", rownames(group_counts)[r[1]], " × ", colnames(group_counts)[r[2]],
             " (n=", group_counts[r[1], r[2]], ")")
    }), collapse = "\n")
  )
}

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n=== Input summary ===")
message("Cells      : ", ncol(obj))
message("Genes      : ", nrow(obj))
message("Cell types : ", paste(sort(unique(obj@meta.data[[opt$cell_type_col]])), collapse = ", "))
message("Samples    : ", paste(sort(unique(obj@meta.data[[opt$sample_col]])),    collapse = ", "))
message("Conditions : ", paste(sort(unique(obj@meta.data[[opt$condition_col]])), collapse = ", "))
message("=====================\n")

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$output_rds)
message("[input_loader] Saved validated object to: ", opt$output_rds)
