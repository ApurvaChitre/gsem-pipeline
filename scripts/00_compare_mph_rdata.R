#!/usr/bin/env Rscript

# Compare GenomicSEM input objects produced by MPH between two runs.
# This is a small debug helper to make sure dimnames/rownames are identical.
#
# Example:
#   Rscript scripts/00_compare_mph_rdata.R \
#     --rdata_a /path/runA/gsem/MPH_genomicSEM.RData \
#     --rdata_b /path/runB/gsem/MPH_genomicSEM.RData

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

rdata_a <- get_arg("--rdata_a")
rdata_b <- get_arg("--rdata_b")
if (is.null(rdata_a) || is.null(rdata_b)) {
  stop("Usage: --rdata_a <file> --rdata_b <file>")
}

load_one <- function(path) {
  objs <- load(path)
  if ("mph_out" %in% objs) return(get("mph_out"))
  get(objs[1])
}

a <- load_one(rdata_a)
b <- load_one(rdata_b)

cmp_dimnames <- function(name, A, B) {
  da <- dimnames(A)
  db <- dimnames(B)
  ok <- identical(da, db)
  list(name=name, identical=ok,
       dimA=paste(dim(A), collapse="x"),
       dimB=paste(dim(B), collapse="x"),
       rownamesA=!is.null(rownames(A)),
       rownamesB=!is.null(rownames(B)))
}

out <- list(
  S = cmp_dimnames("S", a$S, b$S),
  V = cmp_dimnames("V", a$V, b$V)
)

dt <- rbindlist(lapply(out, as.data.table), fill=TRUE)
print(dt)

if (any(dt$identical == FALSE)) {
  message("\nNOT IDENTICAL: At least one of S/V differs in dimnames between the two runs.")
  quit(status = 1)
} else {
  message("\nOK: S and V dimnames are identical.")
}
