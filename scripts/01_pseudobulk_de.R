#!/usr/bin/env Rscript
# 01_pseudobulk_de.R
# Pseudobulk differential expression per cell type using DESeq2 or edgeR.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",         type = "character", help = "Path to validated Seurat .rds"),
  make_option("--method",        type = "character", default = "DESeq2",
              help = "DE method: DESeq2 | edgeR [default: %default]"),
  make_option("--ref_level",     type = "character", default = "",
              help = "Reference condition level (empty = first alphabetically)"),
  make_option("--covariates",    type = "character", default = "",
              help = "Comma-separated covariate column names for the model"),
  make_option("--padj_cutoff",   type = "double",    default = 0.05),
  make_option("--log2fc_cutoff", type = "double",    default = 1.0),
  make_option("--output_dir",    type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$input), !is.null(opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
message("[pseudobulk_de] Loading: ", opt$input)
obj <- readRDS(opt$input)

meta    <- obj@meta.data
ct_col  <- attr(obj, "cell_type_col")  %||% guess_col(meta, c("cell_type","celltype","CellType"))
smp_col <- attr(obj, "sample_col")     %||% guess_col(meta, c("sample_id","sample","Sample"))
cnd_col <- attr(obj, "condition_col")  %||% guess_col(meta, c("condition","Condition","group"))

# Helper: pick first matching column name
`%||%` <- function(a, b) if (!is.null(a) && nchar(a) > 0) a else b

guess_col <- function(meta, candidates) {
  hit <- intersect(candidates, colnames(meta))
  if (length(hit) == 0) stop("Cannot auto-detect column. Pass explicit column names.")
  hit[1]
}

# Re-detect from metadata if attrs not set
if (!exists("ct_col")  || is.null(ct_col))  ct_col  <- guess_col(meta, c("cell_type","celltype","CellType","cluster"))
if (!exists("smp_col") || is.null(smp_col)) smp_col <- guess_col(meta, c("sample_id","sample","Sample"))
if (!exists("cnd_col") || is.null(cnd_col)) cnd_col <- guess_col(meta, c("condition","Condition","group","Group"))

# Set reference level
conditions <- unique(meta[[cnd_col]])
if (nchar(opt$ref_level) == 0) {
  ref_level <- sort(conditions)[1]
  message("[pseudobulk_de] Reference level set to: ", ref_level)
} else {
  ref_level <- opt$ref_level
  if (!ref_level %in% conditions) stop("ref_level '", ref_level, "' not found in condition column.")
}
meta[[cnd_col]] <- relevel(factor(meta[[cnd_col]]), ref = ref_level)

# Covariates
covariates <- character(0)
if (nchar(opt$covariates) > 0) {
  covariates <- trimws(strsplit(opt$covariates, ",")[[1]])
  missing_cov <- setdiff(covariates, colnames(meta))
  if (length(missing_cov) > 0) stop("Covariates not in metadata: ", paste(missing_cov, collapse = ", "))
}

# ── Pseudobulk aggregation ────────────────────────────────────────────────────
message("[pseudobulk_de] Aggregating pseudobulk counts per cell_type × sample...")

cell_types <- sort(unique(meta[[ct_col]]))
all_results <- list()

for (ct in cell_types) {
  message("  Processing: ", ct)
  cells <- rownames(meta)[meta[[ct_col]] == ct]
  sub   <- subset(obj, cells = cells)
  sub_meta <- sub@meta.data

  # Sum counts per sample
  samples  <- sort(unique(sub_meta[[smp_col]]))
  count_mx <- sapply(samples, function(s) {
    s_cells <- rownames(sub_meta)[sub_meta[[smp_col]] == s]
    if (length(s_cells) < 2) return(rep(NA_real_, nrow(sub)))
    Matrix::rowSums(GetAssayData(sub, slot = "counts")[, s_cells, drop = FALSE])
  })
  colnames(count_mx) <- samples

  # Remove samples with too few cells or NA
  valid_samples <- samples[colSums(is.na(count_mx)) == 0]
  count_mx <- count_mx[, valid_samples, drop = FALSE]

  if (ncol(count_mx) < 2) {
    warning("  Skipping ", ct, ": fewer than 2 valid samples after filtering.")
    next
  }

  # Build colData
  col_data <- sub_meta[!duplicated(sub_meta[[smp_col]]), , drop = FALSE]
  rownames(col_data) <- col_data[[smp_col]]
  col_data <- col_data[valid_samples, , drop = FALSE]
  col_data[[cnd_col]] <- relevel(factor(col_data[[cnd_col]]), ref = ref_level)

  # Model formula
  formula_str <- paste("~", paste(c(covariates, cnd_col), collapse = " + "))
  message("    Formula: ", formula_str)

  tryCatch({
    if (opt$method == "DESeq2") {
      res_df <- run_deseq2(count_mx, col_data, formula_str, cnd_col, ref_level)
    } else if (opt$method == "edgeR") {
      res_df <- run_edger(count_mx, col_data, cnd_col, ref_level, covariates)
    } else {
      stop("Unknown method: ", opt$method)
    }
    res_df$cell_type <- ct
    all_results[[ct]] <- res_df
  }, error = function(e) {
    warning("  Error processing ", ct, ": ", conditionMessage(e))
  })
}

# ── DESeq2 ────────────────────────────────────────────────────────────────────
run_deseq2 <- function(count_mx, col_data, formula_str, cnd_col, ref_level) {
  suppressPackageStartupMessages(library(DESeq2))
  # Remove all-zero genes
  count_mx <- count_mx[rowSums(count_mx) > 0, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mx),
    colData   = col_data,
    design    = as.formula(formula_str)
  )
  dds    <- DESeq(dds, quiet = TRUE)
  levels <- resultsNames(dds)
  # Pick the contrast for condition
  contrast_name <- grep(paste0("^", cnd_col, "_"), levels, value = TRUE)
  if (length(contrast_name) == 0) stop("Cannot find condition contrast in DESeq2 results.")

  res_list <- lapply(contrast_name, function(cn) {
    res <- results(dds, name = cn, independentFiltering = TRUE)
    as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      mutate(contrast = cn)
  })
  bind_rows(res_list) %>%
    rename(log2FoldChange = log2FoldChange, padj = padj, pvalue = pvalue,
           baseMean = baseMean, lfcSE = lfcSE, stat = stat)
}

# ── edgeR ─────────────────────────────────────────────────────────────────────
run_edger <- function(count_mx, col_data, cnd_col, ref_level, covariates) {
  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
  })
  count_mx <- count_mx[rowSums(count_mx) > 0, , drop = FALSE]
  y <- DGEList(counts = count_mx)
  y <- calcNormFactors(y)

  # Design matrix
  formula_str <- paste("~", paste(c(covariates, cnd_col), collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = col_data)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  # Test each non-intercept condition coefficient
  cond_cols <- grep(cnd_col, colnames(design), value = TRUE)
  res_list  <- lapply(cond_cols, function(cc) {
    qlf <- glmQLFTest(fit, coef = cc)
    topTags(qlf, n = Inf, sort.by = "none")$table %>%
      rownames_to_column("gene") %>%
      rename(log2FoldChange = logFC, pvalue = PValue, padj = FDR) %>%
      mutate(baseMean = NA_real_, lfcSE = NA_real_, stat = F,
             contrast = cc)
  })
  bind_rows(res_list)
}

# ── Save results ──────────────────────────────────────────────────────────────
if (length(all_results) == 0) stop("No cell types produced DE results.")

combined <- bind_rows(all_results)

# Per-cell-type CSV
for (ct in names(all_results)) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  write.csv(all_results[[ct]],
            file = file.path(opt$output_dir, paste0("de_", ct_safe, ".csv")),
            row.names = FALSE)
}

# Combined CSV
write.csv(combined,
          file = file.path(opt$output_dir, "de_all_celltypes.csv"),
          row.names = FALSE)

# Significant genes summary
sig <- combined %>%
  filter(!is.na(padj), padj < opt$padj_cutoff, abs(log2FoldChange) >= opt$log2fc_cutoff)
write.csv(sig,
          file = file.path(opt$output_dir, "de_significant.csv"),
          row.names = FALSE)

message("[pseudobulk_de] Significant DE genes (padj<", opt$padj_cutoff,
        ", |log2FC|>=", opt$log2fc_cutoff, "): ", nrow(sig))

# ── Volcano plots ─────────────────────────────────────────────────────────────
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in names(all_results)) {
  df <- all_results[[ct]] %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      sig = case_when(
        padj < opt$padj_cutoff & log2FoldChange >=  opt$log2fc_cutoff ~ "Up",
        padj < opt$padj_cutoff & log2FoldChange <= -opt$log2fc_cutoff ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(padj + 1e-300)
    )

  top_genes <- df %>% filter(sig != "NS") %>%
    arrange(padj) %>% slice_head(n = 20)

  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  contrasts <- unique(df$contrast)

  for (ctr in contrasts) {
    df_ctr  <- filter(df, contrast == ctr)
    top_ctr <- filter(top_genes, contrast == ctr)

    p <- ggplot(df_ctr, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
      geom_point(alpha = 0.4, size = 0.8) +
      geom_text_repel(data = top_ctr, aes(label = gene), size = 2.5,
                      max.overlaps = 20, show.legend = FALSE) +
      scale_color_manual(values = c(Up = "#e74c3c", Down = "#3498db", NS = "grey70")) +
      geom_vline(xintercept = c(-opt$log2fc_cutoff, opt$log2fc_cutoff),
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(opt$padj_cutoff),
                 linetype = "dashed", color = "grey50") +
      labs(title = paste0(ct, "\n", ctr),
           x = "log2 Fold Change", y = "-log10(padj)", color = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(size = 10))

    ctr_safe <- gsub("[^A-Za-z0-9_]", "_", ctr)
    ggsave(file.path(plot_dir, paste0("volcano_", ct_safe, "_", ctr_safe, ".png")),
           p, width = 6, height = 5, dpi = 150)
  }
}

message("[pseudobulk_de] Volcano plots saved to: ", plot_dir)
message("[pseudobulk_de] Done.")
