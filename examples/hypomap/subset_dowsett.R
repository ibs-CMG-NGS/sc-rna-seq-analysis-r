#!/usr/bin/env Rscript
# Build a small, comparison-ready subset from the HypoMap atlas.
#
# Keeps one dataset (default: Dowsett10xnuc), removes cells flagged by the
# authors, adds the standard column names expected by the downstream modules
# (cell_type / sample_id / condition), and drops cell types that do not have
# enough cells in every sample.
#
# Usage:
#   Rscript examples/hypomap/subset_dowsett.R \
#     --input  data/hypomap/hypoMap.rds \
#     --output data/hypomap/dowsett_subset.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option("--input",            type = "character",
              default = "data/hypomap/hypoMap.rds"),
  make_option("--output",           type = "character",
              default = "data/hypomap/dowsett_subset.rds"),
  make_option("--dataset",          type = "character", default = "Dowsett10xnuc"),
  make_option("--group_column",     type = "character", default = "C66_named"),
  make_option("--sample_column",    type = "character", default = "Sample_ID"),
  make_option("--condition_column", type = "character", default = "Diet"),
  make_option("--condition_levels", type = "character",
              default = "Normal chow,Fasted"),
  make_option("--exclude_column",   type = "character", default = "Author_Exclude"),
  make_option("--exclude_values",   type = "character", default = "yes"),
  make_option("--min_cells",        type = "integer",   default = 10L,
              help = "Minimum cells per cell type in EVERY sample [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(file.exists(opt$input))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

message("[subset] Loading ", opt$input)
message("[subset] This requires ~25 GB RAM and may take a few minutes.")
obj <- readRDS(opt$input)
if (!inherits(obj, "Seurat")) stop("[subset] input is not a Seurat object")

md <- obj@meta.data
for (col in c("Dataset", opt$group_column, opt$sample_column,
              opt$condition_column)) {
  if (!col %in% colnames(md)) {
    stop("[subset] metadata column not found: ", col)
  }
}

keep <- md$Dataset == opt$dataset
if (nzchar(opt$exclude_column) && opt$exclude_column %in% colnames(md)) {
  excl <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  excl <- excl[nzchar(excl)]
  keep <- keep & (is.na(md[[opt$exclude_column]]) |
                    !(md[[opt$exclude_column]] %in% excl))
}
if (!any(keep)) stop("[subset] no cells left after dataset/exclusion filter")

md <- md[keep, , drop = FALSE]
message("[subset] ", opt$dataset, ": ", nrow(md), " cells after exclusion")

cond_levels <- trimws(strsplit(opt$condition_levels, ",")[[1]])
cond_levels <- cond_levels[nzchar(cond_levels)]
observed <- unique(as.character(md[[opt$condition_column]]))
if (!all(cond_levels %in% observed)) {
  stop("[subset] requested condition levels not found: ",
       paste(setdiff(cond_levels, observed), collapse = ", "),
       "\n  observed: ", paste(observed, collapse = ", "))
}

# Standard columns used by 01_pseudobulk_de.R and 01b_composition.R.
md$cell_type <- as.character(md[[opt$group_column]])
md$sample_id <- as.character(md[[opt$sample_column]])
md$condition <- factor(as.character(md[[opt$condition_column]]),
                       levels = cond_levels)
obj@meta.data <- md

# Keep only cell types with >= min_cells in EVERY sample.
tab <- table(md$cell_type, md$sample_id)
keep_types <- rownames(tab)[rowSums(tab >= opt$min_cells) == ncol(tab)]
if (!length(keep_types)) {
  stop("[subset] no cell type has >= ", opt$min_cells,
       " cells in every sample")
}
message("[subset] cell types kept: ", length(keep_types), " / ", nrow(tab))

cell_ids <- rownames(md)[md$cell_type %in% keep_types]
obj <- subset(obj, cells = cell_ids)
md <- obj@meta.data
md$cell_type <- droplevels(factor(md$cell_type))
md$sample_id <- droplevels(factor(md$sample_id))
md$condition <- factor(md$condition, levels = cond_levels)
obj@meta.data <- md
rm(md, tab)
gc(verbose = FALSE)

cat("\n=== Dowsett subset ===\n")
cat("Cells      :", ncol(obj), "\n")
cat("Genes      :", nrow(obj), "\n")
cat("Cell types :", length(unique(obj$cell_type)), "\n")
cat("Samples    :", length(unique(obj$sample_id)), "\n")
cat("Conditions :\n")
print(table(obj$condition))
cat("Cells per sample:\n")
print(table(obj$sample_id, obj$condition))
cat("======================\n\n")

saveRDS(obj, opt$output)
message("[subset] Saved: ", opt$output, " (",
        round(file.info(opt$output)$size / 1024^2, 1), " MB)")
