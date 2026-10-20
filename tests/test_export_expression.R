#!/usr/bin/env Rscript
# Smoke test for scripts/07_export_expression.R.
# Builds a small synthetic Seurat object, runs the export twice
# (overall and by_dataset), and checks aggregates against an
# independent recomputation. No external fixture files needed.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(arrow)
  library(jsonlite)
})

root <- normalizePath(".")
script <- file.path(root, "scripts", "07_export_expression.R")
stopifnot(file.exists(script))

tmp <- file.path(tempdir(), paste0("exprexport_", as.integer(Sys.time())))
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

set.seed(42)
genes <- sprintf("Gene%03d", 1:200)
n_per <- c(Neuron = 200L, Astro = 200L, Microglia = 200L, Rare = 5L)
cell_type <- rep(names(n_per), n_per)
n_cells <- length(cell_type)
datasets <- ifelse(cell_type == "Rare", "D1", rep(c("D1", "D2"), length.out = n_cells))

counts <- matrix(rbinom(length(genes) * n_cells, size = 1, prob = 0.05),
                 nrow = length(genes), dimnames = list(genes, NULL))
bump <- function(gene, ct, lo, hi) {
  idx <- which(cell_type == ct)
  counts[gene, idx] <<- sample(lo:hi, length(idx), replace = TRUE)
}
bump("Gene001", "Neuron", 3L, 9L)
bump("Gene002", "Astro", 4L, 12L)
bump("Gene003", "Microglia", 2L, 6L)

exclude <- rep("no", n_cells)
exclude[which(cell_type == "Neuron")[1:10]] <- "yes"

obj <- CreateSeuratObject(counts = counts, project = "fixture")
obj <- SetAssayData(obj, assay = "RNA", layer = "data",
                    new.data = log1p(counts))
obj$cell_type <- cell_type
obj$Dataset <- datasets
obj$Author_Exclude <- exclude
fixture <- file.path(tmp, "fixture.rds")
saveRDS(obj, fixture)

run_export <- function(out, extra = character()) {
  log <- file.path(tmp, paste0(basename(out), ".log"))
  args <- c(script,
            "--input", fixture,
            "--output_dir", out,
            "--cell_type_column", "cell_type",
            "--min_cells", "10",
            extra)
  status <- system2("Rscript", args, stdout = log, stderr = log)
  if (status != 0) {
    cat(readLines(log), sep = "\n")
    stop("export failed for ", out)
  }
}

check <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg)
  cat("  ok:", msg, "\n")
}

# --- overall (no dataset split) -------------------------------------------
cat("[overall]\n")
out1 <- file.path(tmp, "overall")
run_export(out1)
long <- read_parquet(file.path(out1, "expression_long.parquet"))
wide <- read_parquet(file.path(out1, "expression_wide_mean.parquet"))
cts  <- read_parquet(file.path(out1, "cell_type_summary.parquet"))
gidx <- read_parquet(file.path(out1, "gene_index.parquet"))
sk   <- read_parquet(file.path(out1, "skipped_groups.parquet"))
mf   <- fromJSON(file.path(out1, "manifest.json"))

kept <- exclude == "no" & cell_type != "Rare"
exp_overall <- sum(vapply(unique(cell_type[kept]), function(ct) {
  idx <- which(kept & cell_type == ct)
  sum(rowSums(counts[, idx, drop = FALSE] > 0) > 0)
}, integer(1)))
check(nrow(long) == exp_overall,
      sprintf("long rows == detected gene-celltype pairs (%d)", exp_overall))
check(nrow(wide) == length(genes) && ncol(wide) == 4L, "wide is 200 x (gene_id + 3 cell types)")
check(nrow(cts) == 3L && all(sort(cts$cell_type) == c("Astro", "Microglia", "Neuron")),
      "cell_type_summary has 3 valid cell types")
check(sum(cts$n_cells) == sum(kept), "summary n_cells == kept cells (590)")
check(nrow(sk) == 1L && sk$group[1] == "Rare" && sk$n_cells[1] == 5L,
      "Rare (5 cells) listed in skipped_groups")
check(nrow(gidx) == length(genes), "gene_index covers all input genes")
check(mf$n_cells_input == 595L && mf$n_cells_used == 590L,
      "manifest n_cells_input=595, n_cells_used=590")
check(identical(unlist(mf$exclusion$values), "yes"), "manifest exclusion values is an array")

for (g in c("Gene001", "Gene002", "Gene003")) {
  ct <- c(Gene001 = "Neuron", Gene002 = "Astro", Gene003 = "Microglia")[[g]]
  idx <- which(cell_type == ct & kept)
  d <- log1p(counts)[g, idx]
  c0 <- counts[g, idx]
  row <- long[long$gene_id == g & long$cell_type == ct, , drop = FALSE]
  check(nrow(row) == 1L, paste(g, ct, "row present"))
  check(isTRUE(all.equal(row$n_cells[1], length(idx))), paste(g, ct, "n_cells"))
  check(isTRUE(all.equal(row$mean_expression[1], mean(d), tolerance = 1e-10)),
        paste(g, ct, "mean_expression matches recomputation"))
  check(isTRUE(all.equal(row$pct_expressing[1], 100 * mean(c0 > 0), tolerance = 1e-10)),
        paste(g, ct, "pct_expressing matches recomputation"))
}

check(!any(is.na(long$mean_expression)) && !any(is.na(long$pct_expressing)),
      "no NA in long metrics")
check(isTRUE(all.equal(long$pct_expressing, 100 * long$n_expressing / long$n_cells,
                       tolerance = 1e-10)), "pct = 100 * n_expressing / n_cells")
top <- gidx[gidx$gene_id == "Gene001", ]
check(top$top_cell_type[1] == "Neuron", "gene_index top_cell_type for Gene001 is Neuron")

# --- by_dataset ------------------------------------------------------------
cat("[by_dataset]\n")
out2 <- file.path(tmp, "by_dataset")
run_export(out2, c("--by_dataset", "TRUE", "--dataset_column", "Dataset"))
long2 <- read_parquet(file.path(out2, "expression_long.parquet"))
cts2  <- read_parquet(file.path(out2, "cell_type_summary.parquet"))
sk2   <- read_parquet(file.path(out2, "skipped_groups.parquet"))

check(nrow(cts2) == 6L, "6 cell_type x dataset groups with >= 10 cells")
check(all(c("dataset", "top_dataset", "top_dataset_frac") %in% colnames(cts2)),
      "by_dataset summary has dataset columns")
check(sum(cts2$n_cells) == sum(kept), "by_dataset summary n_cells == kept cells (590)")
gkey <- paste(cell_type, datasets, sep = " || ")
exp_bd <- sum(vapply(unique(gkey[kept]), function(k) {
  idx <- which(kept & gkey == k)
  sum(rowSums(counts[, idx, drop = FALSE] > 0) > 0)
}, integer(1)))
check(nrow(long2) == exp_bd,
      sprintf("by_dataset long rows == detected gene-group pairs (%d)", exp_bd))
check(nrow(sk2) == 1L && sk2$cell_type[1] == "Rare" && sk2$dataset[1] == "D1",
      "Rare || D1 listed in skipped_groups with cell_type/dataset columns")

idx <- which(cell_type == "Neuron" & datasets == "D2" & kept)
row <- long2[long2$gene_id == "Gene001" & long2$cell_type == "Neuron" &
             long2$dataset == "D2", , drop = FALSE]
check(nrow(row) == 1L, "by_dataset Gene001 / Neuron / D2 row present")
check(isTRUE(all.equal(row$n_cells[1], length(idx))), "by_dataset n_cells matches")
check(isTRUE(all.equal(row$mean_expression[1], mean(log1p(counts)["Gene001", idx]),
                       tolerance = 1e-10)), "by_dataset mean matches recomputation")

# --- viewer bundle ----------------------------------------------------------
cat("[viewer bundle]\n")
bundle <- file.path(tmp, "bundle")
log <- file.path(tmp, "bundle.log")
status <- system2("Rscript",
  c(file.path(root, "scripts", "08_build_viewer_bundle.R"),
    "--export_dir", out1,
    "--bundle_dir", bundle,
    "--prefix", "fixture",
    "--alias", "Fixture_expression",
    "--cell_type_level", "test",
    "--fallback_multi_group", "TRUE"),
  stdout = log, stderr = log)
if (status != 0) {
  cat(readLines(log), sep = "\n")
  stop("bundle build failed")
}
mf_b <- fromJSON(file.path(bundle, "seqviewer_manifest.json"))
check(mf_b$datasets$dataset_type == "expression_summary", "bundle dataset_type")
check(mf_b$datasets$expression_schema$layout == "wide", "bundle layout is wide")
check(mf_b$datasets$expression_schema$pct_scale == "0-100", "bundle declares pct_scale 0-100")
check(mf_b$datasets$file_path == "fixture_wide_mean.parquet",
      "bundle file_path is a bare filename")
cf <- unlist(mf_b$datasets$expression_schema$companion_files)
check(all(file.exists(file.path(bundle, "datasets", cf))),
      "bundle companion files exist")
wb <- read_parquet(file.path(bundle, "datasets", "fixture_wide_mean.parquet"))
check(identical(names(wb)[1:2], c("gene_id", "gene_name")),
      "bundle wide starts with gene_id, gene_name")
mg <- read_parquet(file.path(bundle, "datasets", "fixture_wide_mean_multi_group.parquet"))
check(!"gene_name" %in% names(mg), "multi_group fallback excludes gene_name")
check(file.exists(file.path(bundle, "seqviewer_manifest.multi_group.json")),
      "multi_group fallback manifest exists")

cat("\nPASS test_export_expression\n")
