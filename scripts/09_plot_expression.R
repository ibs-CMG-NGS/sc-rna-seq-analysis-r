#!/usr/bin/env Rscript
# Static summary figures for an expression export (07_export_expression.R).
# Reads the parquet outputs only; does not reload the Seurat object.

suppressPackageStartupMessages({
  library(arrow)
  library(optparse)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(scales)
})

option_list <- list(
  make_option("--input_dir",    type = "character", help = "07_export_expression.R output dir"),
  make_option("--output_dir",   type = "character", default = ""),
  make_option("--top_genes",    type = "integer",   default = 50L),
  make_option("--marker_top_n", type = "integer",   default = 40L),
  make_option("--genes",        type = "character", default = ""),
  make_option("--max_cell_types", type = "integer", default = 0L),
  make_option("--min_pct",      type = "numeric",   default = 10),
  make_option("--min_mean",     type = "numeric",   default = 0.05),
  make_option("--dpi",          type = "integer",   default = 150L),
  make_option("--pdf",          type = "character", default = "TRUE")
)
opt <- parse_args(OptionParser(option_list = option_list))
as_bool <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES", "Y")
make_pdf <- as_bool(opt$pdf)
stopifnot(nzchar(opt$input_dir), dir.exists(opt$input_dir))
if (!nzchar(opt$output_dir)) opt$output_dir <- file.path(opt$input_dir, "plots")
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

read_in <- function(name) {
  p <- file.path(opt$input_dir, name)
  if (!file.exists(p)) stop("[plot] missing ", p)
  arrow::read_parquet(p)
}
wm  <- read_in("expression_wide_mean.parquet")
wp  <- read_in("expression_wide_pct.parquet")
gidx<- read_in("gene_index.parquet")
cs  <- read_in("cell_type_summary.parquet")
mf  <- fromJSON(file.path(opt$input_dir, "manifest.json"))
sk_path <- file.path(opt$input_dir, "skipped_groups.parquet")
sk <- if (file.exists(sk_path)) arrow::read_parquet(sk_path) else NULL
cds_path <- file.path(opt$input_dir, "cell_type_dataset_counts.parquet")
cds <- if (file.exists(cds_path)) arrow::read_parquet(cds_path) else NULL

id_cols <- intersect(c("gene_id", "gene_name"), names(wm))
ct_cols <- setdiff(names(wm), id_cols)
m <- as.matrix(wm[, ct_cols, drop = FALSE]); rownames(m) <- wm$gene_id
p <- as.matrix(wp[, ct_cols, drop = FALSE]); rownames(p) <- wp$gene_id

# keep gene_index aligned with the matrix rows
gidx <- gidx[match(rownames(m), gidx$gene_id), , drop = FALSE]

cs_key <- if ("dataset" %in% names(cs)) {
  paste0(cs$cell_type, " || ", cs$dataset)
} else {
  cs$cell_type
}
ord <- match(cs_key[order(-cs$n_cells)], colnames(m))
ord <- ord[!is.na(ord)]
m <- m[, ord, drop = FALSE]
p <- p[, ord, drop = FALSE]

n <- nrow(m)
top1 <- max.col(m, ties.method = "first")
t1 <- m[cbind(seq_len(n), top1)]
m2 <- m
m2[cbind(seq_len(n), top1)] <- -Inf
t2 <- m2[cbind(seq_len(n), max.col(m2, ties.method = "first"))]
spec <- t1 - t2
specr <- (t1 - t2) / (t1 + 1e-6)
base_marker <- gidx$max_pct_expressing >= opt$min_pct &
               gidx$max_mean_expression >= opt$min_mean &
               gidx$n_cell_types_detected >= 2
base_enrich <- base_marker &
               gidx$n_cell_types_detected <= 0.9 * ncol(m)

short_ct <- function(x) {
  x <- sub("^C[0-9]+-[0-9]+: ", "", x)
  sub(" \\|\\| ", " | ", x)
}
gene_label <- function(ids) {
  x <- gidx$gene_name[match(ids, gidx$gene_id)]
  bad <- is.na(x) | !nzchar(x)
  x[bad] <- ids[bad]
  make.unique(x)
}
hm_breaks <- seq(-2.5, 2.5, length.out = 101)
hm_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
save_heatmap <- function(z, file_base, width, height, cluster_rows,
                         main) {
  pheatmap(z, cluster_rows = cluster_rows, cluster_cols = FALSE,
           labels_col = short_ct(colnames(z)),
           fontsize_row = 7, fontsize_col = 5, angle_col = 90,
           breaks = hm_breaks, color = hm_cols, main = main,
           filename = paste0(file_base, ".png"),
           width = width, height = height, res = opt$dpi)
  if (make_pdf) {
    pheatmap(z, cluster_rows = cluster_rows, cluster_cols = FALSE,
             labels_col = short_ct(colnames(z)),
             fontsize_row = 7, fontsize_col = 5, angle_col = 90,
             breaks = hm_breaks, color = hm_cols, main = main,
             filename = paste0(file_base, ".pdf"),
             width = width, height = height)
  }
}
save_png <- function(p, name, width, height) {
  ggsave(file.path(opt$output_dir, name), p, width = width, height = height,
         dpi = opt$dpi, limitsize = FALSE, bg = "white")
}

# 1. cells per cell type (skipped groups shown in grey) ----------------------
bar_df <- data.frame(cell_type = cs_key, n_cells = cs$n_cells,
                     status = "valid", stringsAsFactors = FALSE)
if (!is.null(sk) && nrow(sk)) {
  sk_key <- if (all(c("cell_type", "dataset") %in% names(sk))) {
    paste0(sk$cell_type, " || ", sk$dataset)
  } else {
    sk$group
  }
  bar_df <- rbind(bar_df, data.frame(cell_type = sk_key, n_cells = sk$n_cells,
                                     status = "skipped", stringsAsFactors = FALSE))
}
bar_df <- bar_df[order(bar_df$n_cells), ]
bar_df$cell_type <- factor(bar_df$cell_type, levels = bar_df$cell_type)
p_bar <- ggplot(bar_df, aes(x = cell_type, y = n_cells, fill = status)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c(valid = "#3B6EA5", skipped = "#999999")) +
  labs(x = NULL, y = "cells (log10)", fill = NULL,
       title = "Cells per cell type",
       subtitle = paste0("valid: ", nrow(cs), " cell types; skipped: ",
                         if (is.null(sk)) 0 else nrow(sk),
                         " | n_cells_used = ", comma(mf$n_cells_used))) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())
save_png(p_bar, "cell_type_n_cells.png",
         width = 9, height = max(5, 0.16 * nrow(bar_df) + 2))

# 2. cell-type-enriched genes (specificity x expression) ---------------------
enrich_top <- integer(0)
if (sum(base_enrich) >= 2 && ncol(m) >= 2) {
  score <- spec * log1p(t1)
  score[!base_enrich] <- -Inf
  enrich_top <- head(order(score, decreasing = TRUE), opt$top_genes)
  hm <- m[enrich_top, , drop = FALSE]
  if (opt$max_cell_types > 0 && ncol(hm) > opt$max_cell_types) {
    hm <- hm[, seq_len(opt$max_cell_types), drop = FALSE]
  }
  z <- t(scale(t(hm)))
  z[!is.finite(z)] <- 0
  z[z > 2.5] <- 2.5; z[z < -2.5] <- -2.5
  rownames(z) <- gene_label(rownames(hm))
  save_heatmap(z, file.path(opt$output_dir, "heatmap_top_enriched_genes"),
               width = max(8, 0.14 * ncol(z) + 5),
               height = max(5, 0.16 * nrow(z) + 3),
               cluster_rows = TRUE,
               main = paste0("Top ", nrow(z), " cell-type-enriched genes (z-score)"))
}

# 3. one marker gene per cell type (diagonal) -------------------------------
marker_tab <- data.frame(gene = rownames(m), ct = colnames(m)[top1],
                         mean = t1, pct = p[cbind(seq_len(n), top1)],
                         specificity = spec, score = spec * log1p(t1),
                         stringsAsFactors = FALSE)
marker_tab <- marker_tab[base_marker, , drop = FALSE]
marker_tab <- marker_tab[order(match(marker_tab$ct, colnames(m)),
                               -marker_tab$score), , drop = FALSE]
marker_tab <- marker_tab[!duplicated(marker_tab$ct), , drop = FALSE]
marker_tab <- marker_tab[!duplicated(marker_tab$gene), , drop = FALSE]

if (nrow(marker_tab) >= 2) {
  rows <- match(marker_tab$gene, rownames(m))
  hm <- m[rows, , drop = FALSE]
  if (opt$max_cell_types > 0 && ncol(hm) > opt$max_cell_types) {
    hm <- hm[, seq_len(opt$max_cell_types), drop = FALSE]
  }
  z <- t(scale(t(hm)))
  z[!is.finite(z)] <- 0
  z[z > 2.5] <- 2.5; z[z < -2.5] <- -2.5
  rownames(z) <- gene_label(rownames(hm))
  save_heatmap(z, file.path(opt$output_dir, "heatmap_marker_genes"),
               width = max(8, 0.14 * ncol(z) + 5),
               height = max(5, 0.16 * nrow(z) + 3),
               cluster_rows = FALSE,
               main = paste0("Top marker gene per cell type (", nrow(z), " genes)"))
  out_mk <- marker_tab
  out_mk$gene_name <- gene_label(out_mk$gene)
  write.csv(out_mk[, c("ct", "gene_name", "gene", "mean", "pct",
                       "specificity", "score")],
            file.path(opt$output_dir, "marker_genes.csv"), row.names = FALSE)
}

# 4. marker dot plot ---------------------------------------------------------
display_n <- if (opt$marker_top_n > 0) min(opt$marker_top_n, ncol(m)) else ncol(m)
display_ct <- colnames(m)[seq_len(display_n)]
genes_req <- trimws(strsplit(opt$genes, ",")[[1]])
genes_req <- genes_req[nzchar(genes_req)]

marker_df <- NULL
if (length(genes_req)) {
  sel_g <- intersect(genes_req, rownames(m))
  if (length(sel_g)) {
    grid <- expand.grid(gene = sel_g, ct = display_ct, stringsAsFactors = FALSE)
    ii <- match(grid$gene, rownames(m)); jj <- match(grid$ct, colnames(m))
    grid$mean <- m[cbind(ii, jj)]
    grid$pct  <- p[cbind(ii, jj)]
    grid$gene_name <- gene_label(grid$gene)
    grid$specificity <- NA_real_
    marker_df <- grid
  } else {
    warning("[plot] requested genes not found: ", paste(genes_req, collapse = ", "))
  }
} else if (nrow(marker_tab)) {
  cand <- marker_tab[marker_tab$ct %in% display_ct, , drop = FALSE]
  cand <- cand[order(match(cand$ct, display_ct)), , drop = FALSE]
  cand <- head(cand, opt$marker_top_n)
  cand$gene_name <- gene_label(cand$gene)
  marker_df <- cand[, c("gene", "gene_name", "ct", "mean", "pct", "specificity")]
}
if (!is.null(marker_df) && nrow(marker_df)) {
  marker_df <- marker_df[marker_df$ct %in% display_ct, , drop = FALSE]
  display_ct <- display_ct[display_ct %in% unique(marker_df$ct)]
  marker_df$label <- ifelse(is.na(marker_df$gene_name) | !nzchar(marker_df$gene_name),
                            marker_df$gene, marker_df$gene_name)
  marker_df$label <- make.unique(marker_df$label)
  marker_df$ct <- factor(marker_df$ct, levels = display_ct)
  marker_df$label <- factor(marker_df$label, levels = rev(unique(marker_df$label)))
  p_dot <- ggplot(marker_df, aes(x = ct, y = label,
                                 size = pct, color = mean)) +
    geom_point() +
    scale_size_continuous(range = c(0.2, 6), limits = c(0, 100)) +
    scale_color_viridis_c(option = "magma") +
    labs(x = NULL, y = NULL, size = "pct (%)", color = "mean",
         title = if (length(genes_req)) "Requested genes" else "Top marker gene per cell type",
         subtitle = paste0(nrow(marker_df), " genes x ", nlevels(marker_df$ct), " cell types")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid.major = element_line(color = "grey92"))
  save_png(p_dot, "top_marker_dotplot.png",
           width = max(8, 0.18 * nlevels(marker_df$ct) + 4),
           height = max(4, 0.20 * nlevels(marker_df$label) + 2.5))
} else {
  warning("[plot] no marker genes passed the filters; dot plot skipped")
}

# 5. gene detection / expression distributions ------------------------------
p1 <- ggplot(gidx, aes(x = n_cell_types_detected)) +
  geom_histogram(binwidth = 1, fill = "#3B6EA5", color = "white") +
  labs(x = "cell types with detection", y = "genes",
       title = "How broadly genes are detected") +
  theme_minimal(base_size = 10)
p2 <- ggplot(gidx[gidx$max_mean_expression > 0, , drop = FALSE],
             aes(x = max_mean_expression)) +
  geom_histogram(bins = 60, fill = "#B2182B", color = "white") +
  scale_x_log10(labels = comma) +
  labs(x = "max mean expression (log10)", y = "genes",
       subtitle = paste0(sum(gidx$max_mean_expression == 0), " genes undetected")) +
  theme_minimal(base_size = 10)
p_det <- p1 | p2
save_png(p_det, "gene_detection_summary.png", width = 11, height = 4.5)

# 6. dataset composition per cell type (by_dataset only) --------------------
if (!is.null(cds) && "dataset" %in% names(cds)) {
  cd <- cds[!is.na(cds$dataset) & nzchar(cds$dataset), , drop = FALSE]
  if (length(unique(cd$dataset)) >= 2) {
    tot <- aggregate(n_cells ~ cell_type, cd, sum)
    keep <- head(tot$cell_type[order(-tot$n_cells)], 30)
    cd <- cd[cd$cell_type %in% keep, , drop = FALSE]
    cd$frac <- cd$n_cells / ave(cd$n_cells, cd$cell_type, FUN = sum)
    cd$cell_type <- factor(cd$cell_type, levels = rev(keep))
    p_ds <- ggplot(cd, aes(x = cell_type, y = frac, fill = dataset)) +
      geom_col(width = 0.8) +
      coord_flip() +
      scale_y_continuous(labels = percent_format()) +
      labs(x = NULL, y = "fraction of cells", fill = NULL,
           title = "Dataset composition per cell type",
           subtitle = paste0(length(unique(cd$dataset)), " datasets; top ", length(keep), " cell types")) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
    save_png(p_ds, "cell_type_dataset_composition.png",
             width = 9, height = max(5, 0.20 * length(keep) + 2))
  }
}

# 7. plot manifest ----------------------------------------------------------
files <- sort(list.files(opt$output_dir, pattern = "\\.(png|pdf|csv)$"))
jsonlite::write_json(list(
  schema_version = "1.0",
  source = normalizePath(opt$input_dir),
  parameters = list(top_genes = opt$top_genes, marker_top_n = opt$marker_top_n,
                    max_cell_types = opt$max_cell_types, min_pct = opt$min_pct,
                    min_mean = opt$min_mean, dpi = opt$dpi),
  features = list(n_enriched_genes = length(enrich_top),
                  n_marker_genes = nrow(marker_tab)),
  files = files
), file.path(opt$output_dir, "plot_manifest.json"),
pretty = TRUE, auto_unbox = TRUE, null = "null")

message("[plot] ", opt$output_dir, " | ", length(files), " files")
