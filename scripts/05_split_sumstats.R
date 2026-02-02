#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  # Supports both:
  #   --flag=value
  #   --flag value
  hit <- grep(paste0('^', flag, '='), args, value = TRUE)
  if (length(hit) > 0) {
    return(sub(paste0('^', flag, '='), '', hit[1]))
  }
  idx <- which(args == flag)
  if (length(idx) > 0) {
    i <- idx[1]
    if (i == length(args)) stop('Missing value for ', flag)
    return(args[i + 1])
  }
  default
}

get_script_dir <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep('^--file=', full, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd()))
  normalizePath(dirname(sub('^--file=', '', file_arg[1])))
}

as_int <- function(x, default) {
  if (is.na(x) || x == '') return(default)
  suppressWarnings({
    v <- as.integer(x)
  })
  if (is.na(v) || v <= 0) stop('Invalid integer for subset_size: ', x)
  v
}

script_dir <- get_script_dir()
base_dir <- get_arg('--base_dir', NA_character_)
if (is.na(base_dir)) {
  # Default: script is in <base_dir>/scripts
  base_dir <- normalizePath(file.path(script_dir, '..'), mustWork = FALSE)
}
base_dir <- normalizePath(base_dir, mustWork = FALSE)

sumstats_rdata <- get_arg('--sumstats_rdata', file.path(base_dir, 'sum_stats_final', 'mySumstatsGSEM.RData'))
out_dir <- get_arg('--out_dir', file.path(base_dir, 'split_sumstats'))
subset_size <- as_int(get_arg('--subset_size', NA_character_), default = 100000)
object_name <- get_arg('--object_name', NA_character_)

if (!file.exists(sumstats_rdata)) {
  stop('sumstats_rdata not found: ', sumstats_rdata)
}

message('Using base_dir: ', base_dir)
message('Reading sumstats_rdata: ', sumstats_rdata)
message('Writing out_dir: ', out_dir)
message('subset_size: ', subset_size)

# Create output directory (YES, it will create it)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load into isolated env so we can safely find the object
E <- new.env(parent = emptyenv())
loaded_names <- load(sumstats_rdata, envir = E)
if (length(loaded_names) == 0) stop('No objects found in: ', sumstats_rdata)

pick_object <- function() {
  if (!is.na(object_name) && object_name != '') {
    if (!exists(object_name, envir = E, inherits = FALSE)) {
      stop('Requested object_name not found in RData: ', object_name)
    }
    return(get(object_name, envir = E, inherits = FALSE))
  }

  if (exists('LocoSumstatsSEM', envir = E, inherits = FALSE)) {
    return(get('LocoSumstatsSEM', envir = E, inherits = FALSE))
  }

  if (length(loaded_names) == 1) {
    return(get(loaded_names[1], envir = E, inherits = FALSE))
  }

  # Otherwise: choose the largest data.frame/data.table object
  candidates <- lapply(loaded_names, function(nm) {
    obj <- get(nm, envir = E, inherits = FALSE)
    if (is.data.frame(obj) || data.table::is.data.table(obj)) {
      return(list(name = nm, n = nrow(obj), obj = obj))
    }
    NULL
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0) {
    stop('Could not find a data.frame/data.table object in RData. Objects: ', paste(loaded_names, collapse = ', '))
  }
  candidates[[which.max(vapply(candidates, `[[`, integer(1), 'n'))]]$obj
}

sumstats <- pick_object()
if (!(is.data.frame(sumstats) || is.data.table(sumstats))) {
  stop('Selected object is not a data.frame/data.table.')
}
setDT(sumstats)

n <- nrow(sumstats)
num_subsets <- ceiling(n / subset_size)

# Optional: write an index file with subset ranges (handy for debugging)
index_dt <- data.table(
  subset_id = seq_len(num_subsets),
  start_row = ((seq_len(num_subsets) - 1) * subset_size) + 1L,
  end_row = pmin(seq_len(num_subsets) * subset_size, n)
)
fwrite(index_dt, file = file.path(out_dir, 'subset_index.tsv'), sep = '\t')

for (i in seq_len(num_subsets)) {
  subset_start <- index_dt$start_row[i]
  subset_end <- index_dt$end_row[i]
  subset <- sumstats[subset_start:subset_end]

  save(subset, file = file.path(out_dir, paste0('sumstats_subset_', i, '.RData')))
}

writeLines(as.character(num_subsets), file.path(out_dir, 'num_SNP_sets.txt'))

message('Done. Wrote ', num_subsets, ' subsets to: ', out_dir)


