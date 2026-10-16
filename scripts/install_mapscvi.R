#!/usr/bin/env Rscript
# install_mapscvi.R
# One-time installation helper for the Phase 0a label-transfer module.
# References: https://github.com/lsteuernagel/mapscvi (MIT)
# NOTE: requires python with scvi-tools >= 0.15 in the active reticulate environment.

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("lsteuernagel/mapscvi", upgrade = "never")
message("[install_mapscvi] mapscvi version: ", as.character(packageVersion("mapscvi")))
