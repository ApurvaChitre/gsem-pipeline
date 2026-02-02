#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

# ----------------------------
# Helpers
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  # Supports both:
  #   --flag=value
  #   --flag value
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) {
    return(sub(paste0("^", flag, "="), "", hit[1]))
  }
  idx <- which(args == flag)
  if (length(idx) > 0) {
    i <- idx[1]
    if (i == length(args)) stop("Missing value for ", flag)
    return(args[i + 1])
  }
  default
}

as_logical <- function(x, default = FALSE) {
  if (is.na(x) || x == "") return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

# Get script dir to support "run-from-anywhere" defaults
get_script_dir <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
}

script_dir <- get_script_dir()

# ----------------------------
# Args (with sensible defaults)
# ----------------------------
base_dir <- get_arg("--base_dir", NA_character_)
if (is.na(base_dir)) {
  # Default: assume script lives under <base_dir>/gsem/scripts or similar
  # Adjust if your layout differs.
  base_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
}
base_dir <- normalizePath(base_dir, mustWork = FALSE)

dict_file     <- get_arg("--dict_file", file.path(base_dir, "pheno_cohort_project_dict.csv"))
sumstats_dir  <- get_arg("--sumstats_dir", file.path(base_dir, "sum_stats_final"))
ref_file      <- get_arg("--ref", file.path(base_dir, "ref_gSEM_frq.txt"))
out_rdata     <- get_arg("--out_rdata", file.path(sumstats_dir, "mySumstatsGSEM.RData"))
pattern       <- get_arg("--pattern", "_formatted.txt$")
allow_partial <- as_logical(get_arg("--allow_partial", "FALSE"), default = FALSE)

# ----------------------------
# Checks
# ----------------------------
if (!dir.exists(base_dir)) stop("base_dir does not exist: ", base_dir)
if (!file.exists(dict_file)) stop("dict_file does not exist: ", dict_file)
if (!dir.exists(sumstats_dir)) stop("sumstats_dir does not exist: ", sumstats_dir)
if (!file.exists(ref_file)) stop("ref file does not exist: ", ref_file)

# ----------------------------
# Load dictionary
# ----------------------------
trait_dict <- fread(dict_file)

needed_cols <- c("pheno_file", "N")
missing_cols <- setdiff(needed_cols, names(trait_dict))
if (length(missing_cols) > 0) {
  stop("Dictionary missing required columns: ", paste(missing_cols, collapse = ", "))
}

# project_name is optional; keep if present (useful for report titles), but do not require it.
if (!("project_name" %in% names(trait_dict))) {
  trait_dict[, project_name := NA_character_]
}

trait_dict <- trait_dict[, .(pheno_file, project_name, N)]

# Normalize trait names to match sumstats filenames
trait_dict[, pheno_file := gsub("^regressedlr_", "", pheno_file)]
trait_dict[, pheno_file := gsub("\\.txt$", "", pheno_file)]
trait_dict[, N := as.numeric(N)]

if (anyNA(trait_dict$N)) {
  stop("Some N values in dictionary are NA or non-numeric. Fix pheno_cohort_project_dict.csv.")
}

# ----------------------------
# Find formatted sumstats files
# ----------------------------
files_full <- list.files(sumstats_dir, pattern = pattern, full.names = TRUE)
if (length(files_full) == 0) {
  stop("No files found in sumstats_dir matching pattern: ", pattern, "\nDir: ", sumstats_dir)
}

dt_files <- data.table(
  file_path  = files_full,
  file_name  = basename(files_full)
)

dt_files[, pheno_file := gsub(pattern, "", file_name)]

# ----------------------------
# Match using dict order (stable / reproducible)
# ----------------------------
matched <- merge(
  trait_dict, dt_files[, .(pheno_file, file_path)],
  by = "pheno_file",
  all.x = TRUE,
  sort = FALSE
)

missing_sumstats <- matched[is.na(file_path), pheno_file]
extra_sumstats   <- setdiff(dt_files$pheno_file, trait_dict$pheno_file)

if (length(missing_sumstats) > 0) {
  msg <- paste0(
    "These traits are in the dictionary but missing sumstats files in sumstats_dir:\n  - ",
    paste(missing_sumstats, collapse = "\n  - ")
  )
  if (!allow_partial) stop(msg) else warning(msg)
}

if (length(extra_sumstats) > 0) {
  warning(
    "These sumstats files exist but are NOT in the dictionary (will be ignored):\n  - ",
    paste(extra_sumstats, collapse = "\n  - ")
  )
}

matched <- matched[!is.na(file_path)]

if (nrow(matched) == 0) stop("After matching, zero traits remain. Check naming and inputs.")

# Final vectors
files        <- matched$file_path
trait.names  <- matched$pheno_file
sample_sizes <- matched$N
se.logit     <- rep(FALSE, length(files))

# ----------------------------
# Run GenomicSEM sumstats
# ----------------------------
message("Using base_dir: ", base_dir)
message("Using dict_file: ", dict_file)
message("Using sumstats_dir: ", sumstats_dir)
message("Using ref: ", ref_file)
message("Traits matched: ", nrow(matched))

LocoSumstatsSEM <- sumstats(
  files       = files,
  ref         = ref_file,
  trait.names = trait.names,
  se.logit    = se.logit,
  maf.filter  = 0
)

# ----------------------------
# Save output
# ----------------------------
out_dir <- dirname(out_rdata)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save(LocoSumstatsSEM, file = out_rdata)
message("Saved: ", out_rdata)

