#!/usr/bin/env Rscript
# 04_cell_communication.R
# Cell-cell communication analysis using CellChat or LIANA.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--input",      type = "character", help = "Path to validated Seurat .rds"),
  make_option("--method",     type = "character", default = "CellChat",
              help = "Method: CellChat | LIANA [default: %default]"),
  make_option("--conditions", type = "character", default = "",
              help = "Comma-separated conditions to compare (empty = all)"),
  make_option("--species",    type = "character", default = "human",
              help = "human | mouse [default: %default]"),
  make_option("--output_dir", type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$input), !is.null(opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
message("[cell_comm] Loading: ", opt$input)
obj  <- readRDS(opt$input)
meta <- obj@meta.data

guess_col <- function(meta, candidates) {
  hit <- intersect(candidates, colnames(meta))
  if (length(hit) == 0) stop("Cannot auto-detect column.")
  hit[1]
}
ct_col  <- guess_col(meta, c("cell_type", "celltype", "CellType", "cluster"))
cnd_col <- guess_col(meta, c("condition", "Condition", "group", "Group"))

# Determine conditions to run
conditions_all <- sort(unique(meta[[cnd_col]]))
if (nchar(opt$conditions) > 0) {
  conditions_use <- trimws(strsplit(opt$conditions, ",")[[1]])
} else {
  conditions_use <- conditions_all
}

message("[cell_comm] Method: ", opt$method)
message("[cell_comm] Conditions: ", paste(conditions_use, collapse = ", "))

# ── Method dispatch ───────────────────────────────────────────────────────────
if (opt$method == "CellChat") {
  run_cellchat(obj, meta, ct_col, cnd_col, conditions_use, opt)
} else if (opt$method == "LIANA") {
  run_liana(obj, meta, ct_col, cnd_col, conditions_use, opt)
} else {
  stop("Unknown method: ", opt$method)
}

# ── CellChat ──────────────────────────────────────────────────────────────────
run_cellchat <- function(obj, meta, ct_col, cnd_col, conditions_use, opt) {
  suppressPackageStartupMessages(library(CellChat))

  db_type <- if (opt$species == "mouse") "mouse" else "human"
  cc_db   <- CellChatDB.human
  if (opt$species == "mouse") cc_db <- CellChatDB.mouse

  cc_objects <- list()

  for (cond in conditions_use) {
    message("[cell_comm] CellChat: ", cond)
    cells <- rownames(meta)[meta[[cnd_col]] == cond]
    sub   <- subset(obj, cells = cells)
    Idents(sub) <- sub@meta.data[[ct_col]]

    cc <- createCellChat(object = sub, group.by = ct_col, assay = "RNA")
    cc@DB <- cc_db

    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- projectData(cc, PPI.human)

    cc <- computeCommunProb(cc)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)

    cc_objects[[cond]] <- cc

    # Save per-condition
    saveRDS(cc, file.path(opt$output_dir, paste0("cellchat_", cond, ".rds")))

    # Interaction count & weight plots
    p1 <- netVisual_circle(cc@net$count,  vertex.weight = table(Idents(sub)),
                           weight.scale = TRUE, label.edge = FALSE,
                           title.name = paste0(cond, " - Number of interactions"))
    p2 <- netVisual_circle(cc@net$weight, vertex.weight = table(Idents(sub)),
                           weight.scale = TRUE, label.edge = FALSE,
                           title.name = paste0(cond, " - Interaction strength"))

    png(file.path(plot_dir, paste0("cellchat_", cond, "_count.png")),
        width = 600, height = 600)
    print(p1); dev.off()
    png(file.path(plot_dir, paste0("cellchat_", cond, "_weight.png")),
        width = 600, height = 600)
    print(p2); dev.off()

    # Pathway heatmap
    p3 <- netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing")
    ggsave(file.path(plot_dir, paste0("cellchat_", cond, "_signaling_heatmap.png")),
           p3, width = 8, height = 6, dpi = 150)

    # Bubble plot (all L-R pairs)
    p4 <- netVisual_bubble(cc, remove.isolate = FALSE) +
      labs(title = paste0(cond, " - L-R interactions"))
    ggsave(file.path(plot_dir, paste0("cellchat_", cond, "_bubble.png")),
           p4, width = 10, height = 8, dpi = 150)

    # Save interaction table
    df_lr <- subsetCommunication(cc)
    write.csv(df_lr, file.path(opt$output_dir, paste0("lr_interactions_", cond, ".csv")),
              row.names = FALSE)
  }

  # Comparative analysis (if 2 conditions)
  if (length(cc_objects) == 2) {
    message("[cell_comm] Running comparative CellChat analysis...")
    cc_merged <- mergeCellChat(cc_objects, add.names = names(cc_objects))
    saveRDS(cc_merged, file.path(opt$output_dir, "cellchat_merged.rds"))

    p_diff <- compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2))
    ggsave(file.path(plot_dir, "cellchat_compare_interactions.png"), p_diff,
           width = 6, height = 5, dpi = 150)
  }

  message("[cell_comm] CellChat complete.")
}

# ── LIANA ─────────────────────────────────────────────────────────────────────
run_liana <- function(obj, meta, ct_col, cnd_col, conditions_use, opt) {
  suppressPackageStartupMessages(library(liana))

  for (cond in conditions_use) {
    message("[cell_comm] LIANA: ", cond)
    cells <- rownames(meta)[meta[[cnd_col]] == cond]
    sub   <- subset(obj, cells = cells)
    Idents(sub) <- sub@meta.data[[ct_col]]

    liana_res <- liana_wrap(sub)
    liana_agg <- liana_aggregate(liana_res)

    write.csv(liana_agg,
              file = file.path(opt$output_dir, paste0("liana_", cond, ".csv")),
              row.names = FALSE)

    # Dot plot
    p <- liana_agg %>%
      filter(aggregate_rank < 0.05) %>%
      liana_dotplot(source_groups = unique(Idents(sub)),
                    target_groups = unique(Idents(sub))) +
      labs(title = paste0(cond, " - LIANA top interactions"))

    ggsave(file.path(plot_dir, paste0("liana_", cond, "_dotplot.png")), p,
           width = 10, height = 8, dpi = 150)
  }

  message("[cell_comm] LIANA complete.")
}

message("[cell_comm] Done. Results in: ", opt$output_dir)
