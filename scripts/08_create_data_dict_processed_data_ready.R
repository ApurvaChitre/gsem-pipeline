#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
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

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n")

# -----------------------------
# Inputs (same mega inputs)
# -----------------------------
processed_in <- get_arg("--processed_in",
                        "/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mega/processed_data_ready.csv"
)

data_dict_in <- get_arg("--data_dict_in",
                        "/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mega/data_dict_mega.csv"
)

# Mega trait column to copy FROM (has regressedlr_ prefix)
source_trait_col <- get_arg("--source_trait_col", "regressedlr_locomotor_mega")

# Old data_dict measure name to replace (NO regressedlr_)
source_dict_measure <- get_arg("--source_dict_measure", "locomotor_mega")

# Remove the old mega trait row from the output data_dict? (default YES)
remove_source_dict_measure <- tolower(get_arg("--remove_source_dict_measure", "true")) %in% c("true","t","1","yes","y")

# -----------------------------
# Outputs (where it writes)
# -----------------------------
processed_out <- get_arg("--processed_out",
                         "/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/metal/metal_report/processed_data_ready.csv"
)

data_dict_out <- get_arg("--data_dict_out",
                         "/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/metal/metal_report/data_dict_mega.csv"
)

# -----------------------------
# Measures to create
#
# data_dict convention: measures are WITHOUT the "regressedlr_" prefix.
# (The processed_data_ready.csv columns typically DO have the "regressedlr_" prefix.)
#
# For convenience, you may pass measures either as:
#   - trait_a,trait_b
#   - regressedlr_trait_a,regressedlr_trait_b
#
# Provide either:
#   --measures=locomotor_metal,locomotor_metal_nogc
# or:
#   --measures=@/path/to/measures.txt   (one per line)
# -----------------------------
read_measures <- function(x) {
  if (is.null(x) || x == "") stop("You must provide --measures=... (comma list) or --measures=@file")
  if (startsWith(x, "@")) {
    f <- sub("^@", "", x)
    if (!file.exists(f)) stop(sprintf("Measures file not found: %s", f))
    v <- readLines(f, warn = FALSE)
  } else {
    v <- strsplit(x, ",", fixed = TRUE)[[1]]
  }
  v <- trimws(v)
  v <- v[v != ""]
  v <- unique(v)
  if (length(v) == 0) stop("No measures found after parsing --measures")
  v
}

measures_raw <- read_measures(get_arg("--measures", NULL))

# Normalize: accept measures with or without a leading regressedlr_ prefix,
# but store measures WITHOUT the prefix (data_dict convention).
measures <- sub("^regressedlr_", "", measures_raw)
measures <- trimws(measures)
measures <- measures[measures != ""]
measures <- unique(measures)

if (length(measures) == 0) stop("No valid measures left after stripping optional regressedlr_ prefix")

msg("DEST: processed_out -> %s", processed_out)
msg("DEST: data_dict_out -> %s", data_dict_out)
msg("Using source_trait_col: %s", source_trait_col)
msg("Measures to create (data_dict; no regressedlr_ prefix): %s", paste(measures, collapse = ", "))
msg("Processed_data_ready columns that will be created: %s", paste(paste0("regressedlr_", measures), collapse = ", "))

# -----------------------------
# 1) processed_data_ready: add regressedlr_<measure> columns copied from mega trait
# -----------------------------
if (!file.exists(processed_in)) stop(sprintf("processed_in not found: %s", processed_in))
processed <- fread(processed_in)

if (!(source_trait_col %in% names(processed))) {
  stop(sprintf(
    "source_trait_col '%s' not found in processed_in. Columns are: %s",
    source_trait_col, paste(names(processed), collapse = ", ")
  ))
}

for (m in measures) {
  new_col <- paste0("regressedlr_", m)
  processed[, (new_col) := get(source_trait_col)]
}

# Drop the original mega trait column in output (keeps only the new regressedlr_<measure> columns):
processed[, (source_trait_col) := NULL]

dir.create(dirname(processed_out), recursive = TRUE, showWarnings = FALSE)
fwrite(processed, processed_out)
msg("Wrote processed_out with %d rows and %d cols", nrow(processed), ncol(processed))

# -----------------------------
# 2) data_dict: create one trait row per measure (NO regressedlr_)
# -----------------------------
if (!file.exists(data_dict_in)) stop(sprintf("data_dict_in not found: %s", data_dict_in))
dd <- fread(data_dict_in, na.strings = c("", "NA"))

# Ensure required columns exist
for (colnm in c("measure", "trait_covariate", "description", "covariates")) {
  if (!(colnm %in% names(dd))) dd[, (colnm) := NA_character_]
}

# Remove any existing rows for the measures we're about to write (avoid dupes)
dd <- dd[!(measure %in% measures)]

# Optionally remove the original mega trait row
if (remove_source_dict_measure) {
  dd <- dd[measure != source_dict_measure]
}

# Add new trait rows
new_rows <- data.table(
  measure = measures,
  trait_covariate = "trait",
  description = "",
  covariates = NA_character_
)

dd <- rbind(dd, new_rows, fill = TRUE)

dir.create(dirname(data_dict_out), recursive = TRUE, showWarnings = FALSE)
fwrite(dd, data_dict_out, na = "NA")
msg("Wrote data_dict_out with %d rows (added %d trait rows)", nrow(dd), nrow(new_rows))

