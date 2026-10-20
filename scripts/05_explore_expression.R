#!/usr/bin/env Rscript
# 05_explore_expression.R
# Exploratory: expression of a gene panel across cell types (no replicates needed).
# Reports n_cells, pct expressing, mean/median/max expression, and rank/specificity.
# Optional per-dataset breakdown to check whether a signal is driven by one study.

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option("--input",        type = "character", help = "Seurat .rds"),
  make_option("--genes",        type = "character", default = "",
              help = "Comma-separated gene symbols"),
  make_option("--group_by",     type = "character", default = "C185_named",
              help = "Metadata column for cell types [default: %default]"),
  make_option("--split_by",     type = "character", default = "Dataset",
              help = "Optional secondary grouping (e.g. Dataset); empty = off"),
  make_option("--min_cells",    type = "integer",   default = 10L,
              help = "Minimum cells per group [default: %default]"),
  make_option("--exclude_col",  type = "character", default = "Author_Exclude",
              help = "Metadata flag column for exclusion; empty = off"),
  make_option("--exclude_vals", type = "character", default = "yes",
              help = "Comma-separated values to exclude [default: %default]"),
  make_option("--assay",        type = "character", default = "RNA"),
  make_option("--top_n",        type = "integer",   default = 20L,
              help = "Top cell types per gene for plots/tables [default: %default]"),
  make_option("--output_dir",   type = "character", default = "output/exploration")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(nzchar(opt$input), file.exists(opt$input))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

get_layer <- function(obj, assay, layer) {
  m <- tryCatch(GetAssayData(obj, assay = assay, layer = layer),
                error = function(e) NULL)
  if (is.null(m)) {
    m <- tryCatch(GetAssayData(obj, assay = assay, slot = layer),
                  error = function(e) NULL)
  }
  if (is.null(m)) stop("[explore] cannot read layer '", layer, "' from assay '", assay, "'")
  m
}

message("[explore] Loading ", opt$input)
obj <- readRDS(opt$input)
if (!inherits(obj, "Seurat")) stop("[explore] input is not a Seurat object")
message("[explore] ", nrow(obj), " genes x ", ncol(obj), " cells")

genes <- trimws(strsplit(opt$genes, ",")[[1]])
genes <- genes[nzchar(genes)]
if (!length(genes)) stop("[explore] --genes is empty")
present <- intersect(genes, rownames(obj))
missing <- setdiff(genes, present)
if (!length(present)) stop("[explore] none of the requested genes are present")
if (length(missing)) message("[explore] WARNING missing genes: ", paste(missing, collapse = ", "))

md <- obj@meta.data
if (!opt$group_by %in% colnames(md)) {
  stop("[explore] group_by column '", opt$group_by, "' not found. Available: ",
       paste(colnames(md), collapse = ", "))
}
if (nzchar(opt$split_by) && !opt$split_by %in% colnames(md)) {
  message("[explore] split_by '", opt$split_by, "' not found; disabling split.")
  opt$split_by <- ""
}

keep <- !is.na(md[[opt$group_by]])
if (nzchar(opt$exclude_col) && opt$exclude_col %in% colnames(md) &&
    nzchar(opt$exclude_vals)) {
  ev <- trimws(strsplit(opt$exclude_vals, ",")[[1]])
  ev <- ev[nzchar(ev)]
  keep <- keep & (is.na(md[[opt$exclude_col]]) | !(md[[opt$exclude_col]] %in% ev))
  message("[explore] exclusion: ", opt$exclude_col, " in {", paste(ev, collapse = ", "),
          "} -> removed ", sum(!keep), " cells")
}
message("[explore] cells retained: ", sum(keep), "/", nrow(md))

data_m   <- get_layer(obj, opt$assay, "data")[present, keep, drop = FALSE]
counts_m <- get_layer(obj, opt$assay, "counts")[present, keep, drop = FALSE]

groups <- as.character(md[[opt$group_by]][keep])
splits <- if (nzchar(opt$split_by)) as.character(md[[opt$split_by]][keep]) else NULL
rm(obj, md); gc(verbose = FALSE)

aggregate_gene <- function(g, key) {
  d <- as.numeric(data_m[g, ])
  c <- as.numeric(counts_m[g, ])
  sp <- split(seq_along(key), key)
  sp <- sp[lengths(sp) > 0]
  data.frame(
    gene           = g,
    group          = names(sp),
    n_cells        = lengths(sp),
    n_expressing   = vapply(sp, function(i) sum(c[i] > 0), integer(1)),
    pct_expressing = vapply(sp, function(i) 100 * mean(c[i] > 0), numeric(1)),
    mean_expr      = vapply(sp, function(i) mean(d[i]), numeric(1)),
    median_expr    = vapply(sp, function(i) stats::median(d[i]), numeric(1)),
    max_expr       = vapply(sp, function(i) max(d[i]), numeric(1)),
    mean_counts    = vapply(sp, function(i) mean(c[i]), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# ── Overall (cell_type) ───────────────────────────────────────────────────────
agg <- do.call(rbind, lapply(present, aggregate_gene, key = groups))
agg <- agg[agg$n_cells >= opt$min_cells, , drop = FALSE]
if (!nrow(agg)) stop("[explore] no group has >= min_cells (", opt$min_cells, ")")

res <- do.call(rbind, lapply(split(agg, agg$gene), function(x) {
  x$rank_by_mean     <- rank(-x$mean_expr, ties.method = "min")
  x$relative_to_max  <- x$mean_expr / max(x$mean_expr)
  x
}))
res <- res[order(res$gene, res$rank_by_mean), ]
write.csv(res, file.path(opt$output_dir, "gene_expression_by_celltype.csv"), row.names = FALSE)

top <- do.call(rbind, lapply(split(res, res$gene), function(x)
  head(x[order(-x$mean_expr), ], opt$top_n)))
write.csv(top, file.path(opt$output_dir, "top_celltypes.csv"), row.names = FALSE)

# ── Per-dataset breakdown (optional) ─────────────────────────────────────────
if (!is.null(splits)) {
  key2 <- paste(groups, splits, sep = " || ")
  agg2 <- do.call(rbind, lapply(present, aggregate_gene, key = key2))
  agg2$cell_type <- sub(" \\|\\| .*$", "", agg2$group)
  agg2$split     <- sub("^.* \\|\\| ", "", agg2$group)
  agg2 <- agg2[agg2$n_cells >= opt$min_cells, , drop = FALSE]
  agg2$group <- NULL
  write.csv(agg2, file.path(opt$output_dir, "gene_expression_by_celltype_dataset.csv"),
            row.names = FALSE)
}

# ── Plots ────────────────────────────────────────────────────────────────────
safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

for (g in unique(top$gene)) {
  d <- top[top$gene == g, ]
  p <- ggplot(d, aes(x = stats::reorder(group, mean_expr), y = mean_expr,
                     fill = pct_expressing)) +
    geom_col(width = 0.8) + coord_flip() +
    scale_fill_gradient(low = "grey85", high = "firebrick3") +
    labs(x = NULL, y = "mean expression (normalized)",
         title = g, subtitle = paste0("top ", nrow(d), " ", opt$group_by),
         fill = "% expressing") +
    theme_bw(base_size = 11)
  ggsave(file.path(plot_dir, paste0("bar_", safe(g), ".png")), p,
         width = 8, height = max(3, 0.22 * nrow(d) + 1.5), dpi = 150)
}

top_ct <- unique(top$group)
ct_max <- tapply(res$mean_expr, res$group, max, na.rm = TRUE)
top_ct <- top_ct[order(-ct_max[top_ct])]
hm <- res[res$group %in% top_ct, ]
m <- tapply(hm$mean_expr, list(hm$gene, hm$group), mean)
z <- t(scale(t(m)))
z[is.na(z)] <- 0
hm_df <- data.frame(
  gene = rep(rownames(z), times = ncol(z)),
  group = rep(colnames(z), each = nrow(z)),
  z = as.numeric(z)
)
p_hm <- ggplot(hm_df, aes(x = group, y = gene, fill = z)) +
  geom_tile() + scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3") +
  labs(x = NULL, y = NULL, title = "Gene expression across top cell types (row z-score)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(plot_dir, "heatmap_top_celltypes.png"), p_hm,
       width = max(8, 0.28 * length(top_ct) + 3), height = max(3, 0.45 * length(present) + 2),
       dpi = 150)

p_dot <- ggplot(res, aes(x = group, y = gene, size = pct_expressing, color = mean_expr)) +
  geom_point() + scale_color_gradient(low = "grey80", high = "firebrick3") +
  labs(x = NULL, y = NULL, size = "% expressing", color = "mean expr",
       title = "Gene panel across cell types") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(plot_dir, "dotplot.png"), p_dot,
       width = max(10, 0.22 * length(unique(res$group)) + 4),
       height = max(3, 0.5 * length(present) + 2), dpi = 150)

# ── Summary ──────────────────────────────────────────────────────────────────
sink(file.path(opt$output_dir, "summary.txt"))
cat("Input:", opt$input, "\n")
cat("Cells used:", sum(keep), "\n")
cat("Group column:", opt$group_by, "(", length(unique(groups)), "groups )\n")
cat("Split column:", if (nzchar(opt$split_by)) opt$split_by else "(none)", "\n")
cat("Genes requested:", paste(genes, collapse = ", "), "\n")
cat("Genes found:", paste(present, collapse = ", "), "\n")
cat("Genes missing:", if (length(missing)) paste(missing, collapse = ", ") else "(none)", "\n")
cat("Min cells per group:", opt$min_cells, "\n\n")
for (g in present) {
  x <- top[top$gene == g, ]
  cat("##", g, "\n")
  print(x[, c("group", "n_cells", "pct_expressing", "mean_expr", "mean_counts")],
        row.names = FALSE)
  cat("\n")
}
sink()

cat("[explore] Done. Outputs in:", opt$output_dir, "\n")
