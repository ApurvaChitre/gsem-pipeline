#!/usr/bin/env Rscript

# If GenomicSEM did not install cleanly via conda, this script attempts
# to install it from GitHub using remotes.
#
# Usage:
#   Rscript env/install_genomicsem_from_github.R

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("remotes package is required. Install via conda (r-remotes) or install.packages('remotes').")
}

message("Installing GenomicSEM from GitHub (GenomicSEM/GenomicSEM)...")
remotes::install_github("GenomicSEM/GenomicSEM", upgrade = "never")

message("Done. Test with: library(GenomicSEM)")
