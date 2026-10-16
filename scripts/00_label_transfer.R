#!/usr/bin/env Rscript
# 00_label_transfer.R
# Phase 0a (optional, v2.2): map cells into a frozen reference model
# (mapscvi / HypoMap) and add predicted labels + provenance.
#
# Scope guard (design decision v2.2):
#   - This module ONLY projects labels from a frozen reference. It does not
#     perform QC, normalization (beyond mapscvi-internal prepare step), or
#     self-clustering annotation.
#   - The downstream cell type column is filled from predictions ONLY when
#     it is absent; existing labels are never overwritten.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",          type = "character",
              help = "Path to input .rds (Seurat) or .h5ad (AnnData)"),
  make_option("--format",         type = "character", default = "seurat",
              help = "Input format: seurat | anndata [default: %default]"),
  make_option("--method",         type = "character", default = "mapscvi",
              help = "Reference mapping method: mapscvi [default: %default]"),
  make_option("--label_column",   type = "character", default = "C66_named",
              help = "Reference metadata column to project (HypoMap: C66_named | C185_named)"),
  make_option("--max_epochs",     type = "integer",   default = 20L,
              help = "scVI mapping epochs [default: %default]"),
  make_option("--suffix",         type = "character", default = "query",
              help = "Internal suffix used by mapscvi [default: %default]"),
  make_option("--cell_type_col",  type = "character", default = "cell_type",
              help = "Downstream pipeline cell type column [default: %default]"),
  make_option("--output_rds",     type = "character", help = "Path to save annotated Seurat .rds"),
  make_option("--provenance_out", type = "character", help = "Path to per-cell provenance CSV")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$input), !is.null(opt$output_rds))

set.seed(1)  # deterministic mapping

# ── AnnData helpers (mirrors input_loader.R) ─────────────────────────────────
py_is_null <- function(x) {
  is.null(x) || inherits(x, "python.builtin.NoneType")
}
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
  vals <- if (inherits(cm, "dgCMatrix")) as.numeric(methods::slot(cm, "x")) else as.numeric(cm)
  if (length(vals) == 0 || any(vals < 0) || max(abs(vals - round(vals))) > 1e-8) {
    stop("[label_transfer] Counts source (", src, ") does not look like raw integer UMI counts.",
         call. = FALSE)
  }
  counts <- Matrix::t(cm)
  rownames(counts) <- var_names
  colnames(counts) <- obs_names
  counts
}

# ── Load query object ────────────────────────────────────────────────────────
message("[label_transfer] Loading query: ", opt$input)
if (opt$format == "seurat") {
  obj <- readRDS(opt$input)
  if (!inherits(obj, "Seurat")) stop("[label_transfer] .rds is not a Seurat object.")
} else if (opt$format == "anndata") {
  suppressPackageStartupMessages(library(reticulate))
  ad     <- reticulate::import("anndata")
  adata  <- ad$read_h5ad(opt$input)
  counts <- read_anndata_counts(adata, "counts")
  obj    <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(adata$obs))
} else {
  stop("[label_transfer] Unknown input format: ", opt$format)
}
message("[label_transfer] Query object: ", ncol(obj), " cells, ", nrow(obj), " genes")

# ── Reference mapping ────────────────────────────────────────────────────────
if (!requireNamespace("mapscvi", quietly = TRUE)) {
  stop("[label_transfer] mapscvi is not installed. Run: Rscript scripts/install_mapscvi.R",
       "\n  (Requires python with scvi-tools >= 0.15 — see environment.yml / README.)",
       call. = FALSE)
}
if (opt$method != "mapscvi") {
  stop("[label_transfer] Only method 'mapscvi' is implemented in Phase 0a.")
}

message("[label_transfer] Mapping into reference model (method=", opt$method,
        ", label_col=", opt$label_column, ", max_epochs=", opt$max_epochs, ")...")
model_path <- system.file("extdata/models/hypoMap_harmonized_scVI_model/", package = "mapscvi")
obj <- mapscvi::map_new_seurat_hypoMap(
  obj,
  suffix     = opt$suffix,
  label_col  = opt$label_column,
  max_epochs = opt$max_epochs
)

# ── Inspect results ──────────────────────────────────────────────────────────
meta <- obj@meta.data
req_cols <- c("predicted", "prediction_probability")
missing_pred <- setdiff(req_cols, colnames(meta))
if (length(missing_pred) > 0) {
  stop("[label_transfer] Mapping did not produce expected columns: ",
       paste(missing_pred, collapse = ", "), call. = FALSE)
}
n_labeled  <- sum(!is.na(meta$predicted))
low_conf_n <- sum(meta$prediction_probability < 0.5, na.rm = TRUE)
message("[label_transfer] Mapped labels: ", n_labeled, "/", ncol(obj),
        " cells ; low confidence (prob<0.5): ", low_conf_n)

# ── Populate downstream cell type column (only if absent) ────────────────────
ct_filled <- FALSE
if (!(opt$cell_type_col %in% colnames(meta))) {
  obj@meta.data[[opt$cell_type_col]] <- obj@meta.data$predicted
  ct_filled <- TRUE
  message("[label_transfer] Filled downstream column '", opt$cell_type_col,
          "' from predictions.")
} else {
  message("[label_transfer] Existing column '", opt$cell_type_col,
          "' kept intact; predictions available as 'predicted'.")
}

# ── Provenance ───────────────────────────────────────────────────────────────
scvi_version <- NA_character_
if (opt$format == "anndata" || requireNamespace("reticulate", quietly = TRUE)) {
  scvi_version <- tryCatch(
    as.character(reticulate::py_to_r(reticulate::import("scvi")$`__version__`)),
    error = function(e) NA_character_
  )
}

per_cell <- data.frame(
  cell_id                = colnames(obj),
  predicted              = obj@meta.data$predicted,
  prediction_probability = obj@meta.data$prediction_probability,
  stringsAsFactors       = FALSE
)

if (!is.null(opt$provenance_out)) {
  dir.create(dirname(opt$provenance_out), recursive = TRUE, showWarnings = FALSE)
  write.csv(per_cell, opt$provenance_out, row.names = FALSE)
  prov_path <- tools::file_path_sans_ext(opt$provenance_out)
  prov <- list(
    method                 = opt$method,
    reference_package      = "mapscvi",
    mapscvi_version        = as.character(packageVersion("mapscvi")),
    seurat_version         = as.character(packageVersion("Seurat")),
    scvi_python_version    = scvi_version,
    model_path             = model_path,
    label_column           = opt$label_column,
    max_epochs             = opt$max_epochs,
    input_path             = opt$input,
    input_format           = opt$format,
    n_cells                = ncol(obj),
    n_labeled              = n_labeled,
    low_confidence_n       = low_conf_n,
    cell_type_col_filled   = ct_filled,
    cell_type_col          = opt$cell_type_col,
    seed                   = 1,
    timestamp              = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  write_yaml_file <- yaml::write_yaml
  write_yaml_file(prov, paste0(prov_path, ".yml"))
  message("[label_transfer] Provenance written to ", paste0(prov_path, ".yml"),
          " and ", opt$provenance_out)
}

# ── Save ─────────────────────────────────────────────────────────────────────
dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$output_rds)
message("[label_transfer] Saved annotated object to: ", opt$output_rds)
message("[label_transfer] Done.")
