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
  make_option("--counts_layer",  type = "character", default = "counts",
              help = "AnnData layer holding raw counts [default: %default]"),
  make_option("--output_rds",    type = "character", help = "Path to save validated Seurat .rds")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$input))      stop("--input is required")
if (is.null(opt$output_rds)) stop("--output_rds is required")

# ── AnnData helpers (v2.0: raw counts priority + sparse matrix support) ───────
py_is_null <- function(x) {
  is.null(x) || inherits(x, "python.builtin.NoneType")
}

# Convert a scipy sparse matrix proxy (reticulate) to Matrix::dgCMatrix.
# scipy stores 0-based indices/indptr — identical convention to Matrix slots.
py_sparse_to_dgc <- function(sp) {
  if (identical(sp$format, "csr")) sp <- sp$tocsc()
  new(
    "dgCMatrix",
    i   = as.integer(reticulate::py_to_r(sp$indices)),
    p   = as.integer(reticulate::py_to_r(sp$indptr)),
    x   = as.numeric(reticulate::py_to_r(sp$data)),
    Dim = as.integer(reticulate::py_to_r(sp$shape))
  )
}

coerce_count_matrix <- function(m) {
  if (inherits(m, c("dgCMatrix", "matrix", "Matrix"))) return(m)
  py_sparse_to_dgc(m)
}

mat_values <- function(m, max_n = 2e6) {
  if (inherits(m, "dgCMatrix")) {
    as.numeric(methods::slot(m, "x"))
  } else {
    as.numeric(m[seq_len(min(max_n, length(m)))])
  }
}

looks_like_counts <- function(vals) {
  length(vals) > 0 && all(vals >= 0) && max(abs(vals - round(vals))) <= 1e-8
}

# Extract gene × cell raw counts from an AnnData object.
# Source priority: explicit layer ("counts") -> adata.raw$X -> adata$X.
# The chosen source must be integer-valued (raw UMI counts) or we stop.
read_anndata_counts <- function(adata, layer_name) {
  var_names <- as.character(reticulate::py_to_r(adata$var_names))
  obs_names <- as.character(reticulate::py_to_r(adata$obs_names))

  src <- "X"
  layer_val <- tryCatch(adata$layers$get(layer_name), error = function(e) NULL)
  if (!py_is_null(layer_val)) {
    cm  <- coerce_count_matrix(layer_val)
    src <- paste0("layers[['", layer_name, "']]")
  } else if (!py_is_null(adata$raw) && !py_is_null(adata$raw$X)) {
    cm  <- coerce_count_matrix(adata$raw$X)
    src <- "raw$X"
  } else {
    cm <- coerce_count_matrix(adata$X)
  }

  vals <- mat_values(cm)
  if (!looks_like_counts(vals)) {
    stop(
      "[input_loader] Selected counts source (", src, ") does not look like raw integer UMI counts.\n",
      "  The matrix appears normalized/log-transformed — pseudobulk aggregation would be invalid.\n",
      "  Re-export the .h5ad with a '", layer_name, "' layer or with raw counts in adata.raw$X.",
      call. = FALSE
    )
  }

  # AnnData stores cells × genes; Seurat expects genes × cells.
  counts <- Matrix::t(cm)
  rownames(counts) <- var_names
  colnames(counts) <- obs_names
  message("[input_loader] Counts source: ", src,
          " (", nrow(counts), " genes x ", ncol(counts), " cells)")
  counts
}

# ── Load ──────────────────────────────────────────────────────────────────────
message("[input_loader] Loading input: ", opt$input)

if (opt$format == "seurat") {
  obj <- readRDS(opt$input)
  if (!inherits(obj, "Seurat")) {
    stop("Object loaded from .rds is not a Seurat object.")
  }

} else if (opt$format == "anndata") {
  suppressPackageStartupMessages(library(reticulate))
  ad    <- reticulate::import("anndata")

  adata  <- ad$read_h5ad(opt$input)
  counts <- read_anndata_counts(adata, opt$counts_layer)

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

# ── Design guard: biological replicates per condition (v2.0) ─────────────────
samples_per_cond <- tapply(
  obj@meta.data[[opt$sample_col]],
  obj@meta.data[[opt$condition_col]],
  function(x) length(unique(x))
)
small_conds <- names(samples_per_cond)[samples_per_cond < 2]
if (length(small_conds) > 0) {
  warning(
    "[input_loader] ", length(small_conds), " condition(s) have < 2 biological samples: ",
    paste(small_conds, collapse = ", "),
    "\n  -> Variance is NOT estimable for these conditions; pseudobulk DE and propeller",
    "\n     results involving them will be unreliable (pseudo-replication risk).",
    "\n  -> They are only acceptable for replicate-free modules (CellChat, trajectory exploration)."
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
