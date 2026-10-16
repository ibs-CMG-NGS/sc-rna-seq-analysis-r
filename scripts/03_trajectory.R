#!/usr/bin/env Rscript
# 03_trajectory.R
# Trajectory / pseudotime analysis using Monocle3 or Slingshot.
# NOTE: Root cell specification requires biological knowledge.
#       root_cell_type and root_cluster should be set deliberately in config.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",          type = "character", help = "Path to validated Seurat .rds"),
  make_option("--method",         type = "character", default = "monocle3",
              help = "Trajectory method: monocle3 | slingshot [default: %default]"),
  make_option("--root_cell_type", type = "character", default = "",
              help = "Cell type to use as trajectory root"),
  make_option("--root_cluster",   type = "character", default = "",
              help = "Cluster within root_cell_type (empty = auto-select centroid)"),
  make_option("--output_dir",     type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$input), !is.null(opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
message("[trajectory] Loading: ", opt$input)
obj <- readRDS(opt$input)
meta <- obj@meta.data

guess_col <- function(meta, candidates) {
  hit <- intersect(candidates, colnames(meta))
  if (length(hit) == 0) stop("Cannot auto-detect column from: ", paste(candidates, collapse = ", "))
  hit[1]
}
ct_col <- guess_col(meta, c("cell_type", "celltype", "CellType", "cluster"))

# ── Validate root ──────────────────────────────────────────────────────────────
if (nchar(opt$root_cell_type) == 0) {
  stop(
    "[trajectory] root_cell_type is required for trajectory analysis.\n",
    "Set trajectory.root_cell_type in your config file.\n",
    "Available cell types: ", paste(sort(unique(meta[[ct_col]])), collapse = ", ")
  )
}
if (!opt$root_cell_type %in% meta[[ct_col]]) {
  stop("[trajectory] root_cell_type '", opt$root_cell_type, "' not found in data.\n",
       "Available: ", paste(sort(unique(meta[[ct_col]])), collapse = ", "))
}
message("[trajectory] Root cell type: ", opt$root_cell_type)

# ── Method dispatch ───────────────────────────────────────────────────────────
if (opt$method == "monocle3") {
  run_monocle3(obj, meta, ct_col, opt)
} else if (opt$method == "slingshot") {
  run_slingshot(obj, meta, ct_col, opt)
} else {
  stop("Unknown trajectory method: ", opt$method)
}

# ── Monocle3 ──────────────────────────────────────────────────────────────────
run_monocle3 <- function(obj, meta, ct_col, opt) {
  suppressPackageStartupMessages({
    library(monocle3)
    library(SeuratWrappers)
  })

  message("[trajectory] Converting to CellDataSet (Monocle3)...")
  cds <- as.cell_data_set(obj)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)

  # Identify root cells: cells in root_cell_type closest to UMAP centroid
  root_cells <- identify_root_cells(cds, obj, ct_col, opt$root_cell_type, opt$root_cluster)
  message("[trajectory] Root cells identified: ", length(root_cells))

  cds <- order_cells(cds, root_cells = root_cells)

  # Save pseudotime
  pt <- pseudotime(cds)
  pt_df <- data.frame(cell = names(pt), pseudotime = pt,
                      cell_type = meta[names(pt), ct_col])
  write.csv(pt_df, file.path(opt$output_dir, "pseudotime_monocle3.csv"), row.names = FALSE)

  # UMAP colored by pseudotime
  p <- plot_cells(cds,
                  color_cells_by     = "pseudotime",
                  label_cell_groups  = FALSE,
                  label_leaves       = FALSE,
                  label_branch_points = FALSE,
                  graph_label_size   = 3) +
    labs(title = paste0("Pseudotime (Monocle3)\nRoot: ", opt$root_cell_type))

  ggsave(file.path(opt$output_dir, "plots", "monocle3_pseudotime.png"), p,
         width = 7, height = 6, dpi = 150)

  # Save CDS
  saveRDS(cds, file.path(opt$output_dir, "monocle3_cds.rds"))
  message("[trajectory] Monocle3 complete.")
}

# ── Slingshot ─────────────────────────────────────────────────────────────────
run_slingshot <- function(obj, meta, ct_col, opt) {
  suppressPackageStartupMessages({
    library(slingshot)
    library(SingleCellExperiment)
  })

  message("[trajectory] Running Slingshot...")

  # Use UMAP embedding
  if (!"umap" %in% names(obj@reductions)) {
    message("[trajectory] No UMAP found; running UMAP first...")
    obj <- RunUMAP(obj, dims = 1:30)
  }

  umap_emb  <- Embeddings(obj, "umap")
  cl_labels <- meta[[ct_col]]
  names(cl_labels) <- rownames(meta)

  start_clus <- opt$root_cell_type

  sce <- SingleCellExperiment(
    assays    = list(counts = GetAssayData(obj, layer = "counts")),
    reducedDims = list(UMAP = umap_emb)
  )

  sds <- slingshot(sce,
                   clusterLabels = cl_labels,
                   reducedDim    = "UMAP",
                   start.clus    = start_clus)

  # Extract pseudotime (average across lineages)
  pt_mat <- slingPseudotime(sds)
  pt_avg <- rowMeans(pt_mat, na.rm = TRUE)

  pt_df <- data.frame(
    cell        = rownames(umap_emb),
    pseudotime  = pt_avg,
    cell_type   = cl_labels[rownames(umap_emb)]
  )
  for (i in seq_len(ncol(pt_mat))) {
    pt_df[[paste0("lineage_", i)]] <- pt_mat[, i]
  }
  write.csv(pt_df, file.path(opt$output_dir, "pseudotime_slingshot.csv"), row.names = FALSE)

  # Plot
  umap_df <- as.data.frame(umap_emb)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$pseudotime <- pt_avg
  umap_df$cell_type  <- cl_labels[rownames(umap_df)]

  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_viridis_c(option = "magma") +
    labs(title = paste0("Pseudotime (Slingshot)\nRoot: ", opt$root_cell_type)) +
    theme_bw(base_size = 11)

  ggsave(file.path(opt$output_dir, "plots", "slingshot_pseudotime.png"), p,
         width = 7, height = 6, dpi = 150)

  saveRDS(sds, file.path(opt$output_dir, "slingshot_sds.rds"))
  message("[trajectory] Slingshot complete.")
}

# ── Root cell helper (Monocle3) ────────────────────────────────────────────────
identify_root_cells <- function(cds, seurat_obj, ct_col, root_cell_type, root_cluster) {
  # Get cells of the root cell type
  root_cell_barcodes <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[ct_col]] == root_cell_type]

  if (nchar(root_cluster) > 0) {
    # Try to narrow to a specific cluster
    cluster_col <- tryCatch(
      { guess_col(seurat_obj@meta.data, c("seurat_clusters", "cluster", "Cluster")); },
      error = function(e) NULL
    )
    if (!is.null(cluster_col)) {
      in_cluster <- rownames(seurat_obj@meta.data)[
        seurat_obj@meta.data[[ct_col]] == root_cell_type &
        as.character(seurat_obj@meta.data[[cluster_col]]) == root_cluster
      ]
      if (length(in_cluster) > 0) root_cell_barcodes <- in_cluster
    }
  }

  # Return the cell closest to the centroid in UMAP
  if ("umap" %in% names(seurat_obj@reductions)) {
    umap <- Embeddings(seurat_obj, "umap")[root_cell_barcodes, , drop = FALSE]
    centroid  <- colMeans(umap)
    distances <- sqrt(rowSums((umap - centroid)^2))
    n_root    <- max(1, round(0.01 * length(root_cell_barcodes)))
    return(names(sort(distances)[seq_len(n_root)]))
  }
  root_cell_barcodes[seq_len(min(5, length(root_cell_barcodes)))]
}

message("[trajectory] Done.")
