#!/usr/bin/env Rscript
# 01b_composition.R
# Cell type compositional analysis using propeller (speckle).

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",      type = "character", help = "Path to validated Seurat .rds"),
  make_option("--output_dir", type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$input), !is.null(opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
message("[composition] Loading: ", opt$input)
obj <- readRDS(opt$input)
meta <- obj@meta.data

# Auto-detect columns
guess_col <- function(meta, candidates) {
  hit <- intersect(candidates, colnames(meta))
  if (length(hit) == 0) stop("Cannot auto-detect column from: ", paste(candidates, collapse = ", "))
  hit[1]
}
ct_col  <- guess_col(meta, c("cell_type", "celltype", "CellType", "cluster"))
smp_col <- guess_col(meta, c("sample_id", "sample", "Sample"))
cnd_col <- guess_col(meta, c("condition", "Condition", "group", "Group"))

message("[composition] Using columns: cell_type=", ct_col, ", sample=", smp_col, ", condition=", cnd_col)

# ── Frequency table ───────────────────────────────────────────────────────────
freq_df <- meta %>%
  group_by(.data[[smp_col]], .data[[ct_col]], .data[[cnd_col]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[smp_col]]) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

colnames(freq_df)[1:3] <- c("sample", "cell_type", "condition")

write.csv(freq_df,
          file = file.path(opt$output_dir, "cell_type_proportions.csv"),
          row.names = FALSE)

# ── Propeller test ────────────────────────────────────────────────────────────
tryCatch({
  suppressPackageStartupMessages(library(speckle))

  message("[composition] Running propeller...")

  # propeller expects: clusters, sample, group vectors (all same length)
  prop_res <- propeller(
    clusters = meta[[ct_col]],
    sample   = meta[[smp_col]],
    group    = meta[[cnd_col]]
  )

  prop_res <- prop_res %>% rownames_to_column("cell_type")
  write.csv(prop_res,
            file = file.path(opt$output_dir, "propeller_results.csv"),
            row.names = FALSE)

  sig_types <- prop_res %>% filter(FDR < 0.05)
  message("[composition] Significant cell types (FDR < 0.05): ", nrow(sig_types))

}, error = function(e) {
  message("[composition] propeller failed: ", conditionMessage(e))
  message("[composition] Falling back to proportions-only output.")
})

# ── Plots ─────────────────────────────────────────────────────────────────────

# Bar plot: stacked proportions per sample
p_bar <- ggplot(freq_df, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(position = "stack") +
  facet_wrap(~ condition, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion", fill = "Cell type",
       title = "Cell type composition per sample") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "composition_bar.png"), p_bar,
       width = max(8, length(unique(freq_df$sample)) * 0.8 + 2),
       height = 5, dpi = 150)

# Box/dot plot: proportion per cell type across conditions
p_box <- ggplot(freq_df, aes(x = condition, y = prop, color = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion", title = "Cell type proportions by condition") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

n_types <- length(unique(freq_df$cell_type))
ggsave(file.path(plot_dir, "composition_by_celltype.png"), p_box,
       width = min(14, ceiling(sqrt(n_types)) * 3),
       height = min(12, ceiling(n_types / ceiling(sqrt(n_types))) * 3),
       dpi = 150)

message("[composition] Plots saved to: ", plot_dir)
message("[composition] Done.")
