#!/usr/bin/env Rscript

# Run GenomicSEM::userGWAS() on ONE split SNP set ("chunk").
# Intended to be called from a Slurm array.

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

mph_rdata   <- get_arg("--mph_rdata")
subset_rdata <- get_arg("--subset_rdata")
out_rdata    <- get_arg("--out_rdata")
dict_file    <- get_arg("--dict_file", default = "")
model_file   <- get_arg("--model_file", default = "")
cores_str    <- get_arg("--cores", default = "1")
parallel_str <- get_arg("--parallel", default = "1")

if (is.null(mph_rdata) || is.null(subset_rdata) || is.null(out_rdata)) {
  stop("Required args: --mph_rdata --subset_rdata --out_rdata")
}

cores <- suppressWarnings(as.integer(cores_str))
if (!is.finite(cores) || cores < 1) cores <- 1
parallel <- parallel_str %in% c("1","true","TRUE","yes","YES")

load_first_object <- function(rdata_path) {
  objs <- load(rdata_path)
  if (length(objs) < 1) stop("No objects found in RData: ", rdata_path)
  get(objs[1])
}

mph_out <- load_first_object(mph_rdata)
sub <- load_first_object(subset_rdata)

# ---- Trait names (prefer mph_out dimnames; fallback to dict) ----
clean_trait_from_pheno_file <- function(x) {
  # Supports MPH dict pheno_file values like:
  #   regressedlr_<trait>.txt
  # and also tolerates paths.
  x <- as.character(x)
  x <- trimws(x)
  x <- basename(x)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  x <- sub("^regressedlr_", "", x)
  x
}

infer_traits_from_covstruc <- function(covstruc) {
  tr <- NULL
  if (!is.null(covstruc$S)) {
    rn <- rownames(covstruc$S)
    cn <- colnames(covstruc$S)
    if (!is.null(rn) && length(rn) > 0 && all(nzchar(rn))) tr <- rn
    if (is.null(tr) && !is.null(cn) && length(cn) > 0 && all(nzchar(cn))) tr <- cn
  }
  if (!is.null(tr)) {
    tr <- tr[!is.na(tr) & tr != ""]
    if (length(tr) == 0) tr <- NULL
  }
  tr
}

traits <- infer_traits_from_covstruc(mph_out)

if (is.null(traits) && nzchar(dict_file) && file.exists(dict_file)) {
  dt <- fread(dict_file)
  if ("trait_name" %in% names(dt)) {
    traits <- as.character(dt$trait_name)
  } else {
    if (!"pheno_file" %in% names(dt)) stop("dict_file must have column 'pheno_file' (or 'trait_name')")
    traits <- clean_trait_from_pheno_file(dt$pheno_file)
  }
}

if (is.null(traits) || length(traits) < 2) {
  stop(
    "Could not determine trait names.\n",
    "Tried mph_out$S dimnames; fallback: --dict_file with column pheno_file (or trait_name).\n",
    "Fix: ensure MPH_genomicSEM.RData has dimnames OR pass --dict_file=/path/to/pheno_cohort_project_dict.csv"
  )
}

# Ensure mph_out$S has dimnames that match the model variables.
if (!is.null(mph_out$S)) {
  mS <- as.matrix(mph_out$S)
  if (nrow(mS) != length(traits) || ncol(mS) != length(traits)) {
    stop(
      "Trait count mismatch:\n",
      "  n_traits (from dict/rownames) = ", length(traits), "\n",
      "  dim(mph_out$S) = ", nrow(mS), " x ", ncol(mS), "\n",
      "Check that your pheno_cohort_project_dict.csv matches the MPH_genomicSEM.RData order and contents."
    )
  }
  dimnames(mS) <- list(traits, traits)
  mph_out$S <- mS
}

if (is.null(mph_out$trait.names)) {
  mph_out$trait.names <- traits
}

# Build model (1-factor userGWAS)
model <- NULL
if (nzchar(model_file) && file.exists(model_file)) {
  model <- paste(readLines(model_file, warn = FALSE), collapse = "\n")
} else {
  model <- paste0("F1 =~ ", paste(traits, collapse = " + "))
}

# Ensure SNP regression is present
if (!grepl("~\\s*SNP", model)) {
  model <- paste(model, "F1 ~ SNP", sep = "\n")
}

# Ensure output directory exists
out_dir <- dirname(out_rdata)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Run userGWAS with clean warning capture
warns <- character(0)
res <- tryCatch({
  withCallingHandlers(
    expr = {
      userGWAS(
        covstruc   = mph_out,
        SNPs       = sub,
        model      = model,
        estimation = "DWLS",
        parallel   = parallel,
        cores      = cores
      )
    },
    warning = function(w) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(w))
      warns <<- c(warns, msg)
      invokeRestart("muffleWarning")
    }
  )
}, error = function(e) {
  # Write a simple sidecar error file
  err_file <- paste0(out_rdata, ".ERROR.txt")
  writeLines(c("userGWAS failed:", conditionMessage(e)), err_file)
  stop(e)
})

warns <- unique(warns)
warn_file <- paste0(out_rdata, ".warnings.txt")
if (length(warns) > 0) {
  writeLines(warns, warn_file)
} else {
  writeLines("(none)", warn_file)
}

save(res, model, traits, warns, file = out_rdata)

message("Done: saved ", out_rdata)
