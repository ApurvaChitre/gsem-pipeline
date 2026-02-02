#!/usr/bin/env Rscript

# Tiny sanity test to confirm warnings are captured and written cleanly.
#
# Example:
#   Rscript scripts/test_warning_capture.R --out_dir /tmp/warn_test

suppressPackageStartupMessages({
  library(data.table)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

out_dir <- get_arg("--out_dir", default = "warn_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

warns <- character(0)
withCallingHandlers(
  expr = {
    warning("This is warning #1, with, commas")
    warning("This is warning #2 with \"quotes\"")
  },
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

warns <- unique(warns)
writeLines(warns, file.path(out_dir, "warnings.txt"))

fwrite(data.table(n_warnings = length(warns)), file.path(out_dir, "summary.tsv"), sep = "\t")

message("Wrote:")
message("  ", file.path(out_dir, "warnings.txt"))
message("  ", file.path(out_dir, "summary.tsv"))
