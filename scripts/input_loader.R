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
# NOTE: reticulate auto-converts scipy sparse matrices to R Matrix objects, and
# that converter is unreliable for non-double dtypes (e.g. int32 UMI counts) —
# it errors with "'x' slot is not of type double". So we never let a sparse
# object cross into R: source selection + CSR extraction happen inside Python
# and only plain numpy arrays (shape/indptr/indices/data) are returned.
ANNDATA_HELPER <- "
import numpy as _np
import scipy.sparse as _sp

def scran_extract_counts(adata, layer_name):
    if layer_name in adata.layers:
        m = adata.layers[layer_name]
        src = \"layers['%s']\" % layer_name
    elif adata.raw is not None and adata.raw.X is not None:
        m = adata.raw.X
        src = 'raw$X'
    else:
        m = adata.X
        src = 'X'
    if hasattr(m, 'tocsr'):
        m = m.tocsr()
    else:
        m = _sp.csr_matrix(_np.asarray(m))
    return {
        'src': src,
        'shape': _np.asarray(m.shape, dtype=_np.int64),      # (cells, genes)
        'indptr': _np.asarray(m.indptr, dtype=_np.int64),
        'indices': _np.asarray(m.indices, dtype=_np.int64),
        'data': _np.asarray(m.data, dtype=_np.float64),
    }
"

looks_like_counts <- function(vals) {
  length(vals) > 0 && all(vals >= 0) && max(abs(vals - round(vals))) <= 1e-8
}

# Extract gene × cell raw counts from an AnnData object.
# Source priority: explicit layer ("counts") -> adata.raw$X -> adata$X.
# The CSR of a cells×genes matrix has the exact storage layout of a
# genes×cells dgCMatrix (indptr = cell pointers, indices = gene indices), so
# Dim is flipped in place — no explicit transpose is needed or used.
read_anndata_counts <- function(adata, layer_name) {
  var_names <- as.character(adata$var_names$to_numpy())
  obs_names <- as.character(adata$obs_names$to_numpy())

  parts  <- reticulate::py$scran_extract_counts(adata, layer_name)
  shape  <- as.integer(parts$shape)                 # c(cells, genes)

  counts <- new(
    "dgCMatrix",
    i   = as.integer(parts$indices),
    p   = as.integer(parts$indptr),
    x   = as.numeric(parts$data),
    Dim = shape[c(2, 1)]                            # genes × cells
  )
  rownames(counts) <- var_names
  colnames(counts) <- obs_names

  if (!looks_like_counts(as.numeric(methods::slot(counts, "x")))) {
    stop(
      "[input_loader] Selected counts source (", parts$src, ") does not look like raw integer UMI counts.\n",
      "  The matrix appears normalized/log-transformed — pseudobulk aggregation would be invalid.\n",
      "  Re-export the .h5ad with a '", layer_name, "' layer or with raw counts in adata.raw$X.",
      call. = FALSE
    )
  }

  message("[input_loader] Counts source: ", parts$src,
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
  reticulate::py_run_string(ANNDATA_HELPER)
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

for (col in c(opt$cell_type_col, opt$sample_col, opt$condition_col)) {
  if (is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- droplevels(obj@meta.data[[col]])
  }
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
