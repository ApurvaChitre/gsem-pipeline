#!/usr/bin/env Rscript

# Fit (and sanity-check) a 1-factor GenomicSEM usermodel using MPH-derived covstruc.
#
# This script is meant to:
#   - confirm the model syntax matches the trait names
#   - surface warnings (Heywood / non-positive definite / etc.) in a clean text file
#   - write a small TSV summary of key fit stats
#
# It does NOT ship multiple alternative models; edit MODEL_FILE if you want to test.

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

mph_rdata <- get_arg("--mph_rdata")
dict_file <- get_arg("--dict_file", "")
out_dir <- get_arg("--out_dir")
model_file <- get_arg("--model_file", "")
estimation <- get_arg("--estimation", "DWLS")
std_lv <- as.logical(get_arg("--std_lv", "TRUE"))

# Preflight / guardrails
# NOTE: This script is used both for:
#   - Step 0 preflight gate (strict-by-config)
#   - Phase A model fit (usually inspect + proceed)

# 1) Strict PD failure (rarely recommended; prefer smooth_diff threshold gate)
fail_on_nonpd <- as.logical(get_arg("--fail_on_nonpd", "FALSE"))

# 2) PD numerical tolerance
pd_tol <- suppressWarnings(as.numeric(get_arg("--pd_tol", "-1e-8")))
if (is.na(pd_tol)) pd_tol <- -1e-8

# 3) Smoothing gate: if S is non-PD, compute nearPD smoothing and measure max abs change.
#    If the max abs change exceeds this threshold, treat as fatal (default 0.025).
fail_on_smooth_exceed <- as.logical(get_arg("--fail_on_smooth_exceed", "TRUE"))
smooth_diff_thresh <- suppressWarnings(as.numeric(get_arg("--smooth_diff_thresh", "0.025")))
if (is.na(smooth_diff_thresh)) smooth_diff_thresh <- 0.025

# 4) Model pathologies
fail_on_heywood <- as.logical(get_arg("--fail_on_heywood", "TRUE"))
fail_on_nonconverged <- as.logical(get_arg("--fail_on_nonconverged", "TRUE"))

if (is.null(mph_rdata) || is.null(out_dir)) {
  stop("Usage: --mph_rdata <MPH_genomicSEM.RData> --out_dir <dir> [--dict_file ...] [--model_file ...]")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Collect fatal issues while still writing as many helpful outputs as possible.
fatal_msgs <- character(0)

load_mph_out <- function(path) {
  objs <- load(path)
  if ("mph_out" %in% objs) return(get("mph_out"))
  get(objs[1])
}

mph_out <- load_mph_out(mph_rdata)

# ---- Check positive-definiteness of the input matrices (S and V) ----
pd_report_one <- function(mat, name, tol = -1e-8) {
  if (is.null(mat)) {
    return(data.table(
      matrix = name,
      n = NA_integer_,
      symmetric = NA,
      na_count = NA_integer_,
      min_eig = NA_real_,
      max_eig = NA_real_,
      n_eig_le_0 = NA_integer_,
      is_pd = NA,
      eigen_error = "(missing)"
    ))
  }

  m <- tryCatch(as.matrix(mat), error = function(e) NULL)
  if (is.null(m)) {
    return(data.table(
      matrix = name,
      n = NA_integer_,
      symmetric = NA,
      na_count = NA_integer_,
      min_eig = NA_real_,
      max_eig = NA_real_,
      n_eig_le_0 = NA_integer_,
      is_pd = NA,
      eigen_error = "could not coerce to matrix"
    ))
  }

  # Force numeric (defensive)
  storage.mode(m) <- "double"

  n <- nrow(m)
  na_count <- sum(is.na(m))
  symm <- tryCatch(isSymmetric(m, tol = 1e-10), error = function(e) NA)

  eig <- tryCatch(
    eigen(m, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) e
  )

  if (inherits(eig, "error")) {
    return(data.table(
      matrix = name,
      n = n,
      symmetric = symm,
      na_count = na_count,
      min_eig = NA_real_,
      max_eig = NA_real_,
      n_eig_le_0 = NA_integer_,
      is_pd = NA,
      eigen_error = conditionMessage(eig)
    ))
  }

  min_eig <- suppressWarnings(min(eig, na.rm = TRUE))
  max_eig <- suppressWarnings(max(eig, na.rm = TRUE))
  n_le_0 <- sum(eig <= 0, na.rm = TRUE)
  is_pd <- is.finite(min_eig) && (min_eig > tol)

  data.table(
    matrix = name,
    n = n,
    symmetric = symm,
    na_count = na_count,
    min_eig = min_eig,
    max_eig = max_eig,
    n_eig_le_0 = n_le_0,
    is_pd = is_pd,
    eigen_error = ""
  )
}

pd_dt <- rbind(
  pd_report_one(mph_out$S, "S", tol = pd_tol),
  pd_report_one(mph_out$V, "V", tol = pd_tol),
  fill = TRUE
)

fwrite(pd_dt, file.path(out_dir, "matrix_pd_check.tsv"), sep = "\t")

nonpd_S <- isTRUE(pd_dt[matrix == "S", is_pd] == FALSE)

# ---- Smoothing check (GenomicSEM convention) ----
# If S is not PD, GenomicSEM may internally apply a nearPD-style smoothing.
# Many workflows proceed if the change from smoothing is small.
smooth_report_one <- function(mat, name, tol = -1e-8, thresh = 0.025) {
  if (is.null(mat)) {
    return(data.table(
      matrix = name,
      smooth_applied = NA,
      method = "Matrix::nearPD",
      max_abs_diff = NA_real_,
      diff_thresh = thresh,
      exceeds_thresh = NA,
      min_eig_original = NA_real_,
      min_eig_smoothed = NA_real_,
      smooth_error = "(missing)"
    ))
  }

  m <- tryCatch(as.matrix(mat), error = function(e) NULL)
  if (is.null(m)) {
    return(data.table(
      matrix = name,
      smooth_applied = NA,
      method = "Matrix::nearPD",
      max_abs_diff = NA_real_,
      diff_thresh = thresh,
      exceeds_thresh = NA,
      min_eig_original = NA_real_,
      min_eig_smoothed = NA_real_,
      smooth_error = "could not coerce to matrix"
    ))
  }

  storage.mode(m) <- "double"
  if (anyNA(m)) {
    return(data.table(
      matrix = name,
      smooth_applied = NA,
      method = "Matrix::nearPD",
      max_abs_diff = NA_real_,
      diff_thresh = thresh,
      exceeds_thresh = NA,
      min_eig_original = NA_real_,
      min_eig_smoothed = NA_real_,
      smooth_error = "matrix contains NA"
    ))
  }

  # Force symmetry to avoid numerical artifacts
  m <- (m + t(m)) / 2

  eig0 <- tryCatch(eigen(m, symmetric = TRUE, only.values = TRUE)$values, error = function(e) e)
  if (inherits(eig0, "error")) {
    return(data.table(
      matrix = name,
      smooth_applied = NA,
      method = "Matrix::nearPD",
      max_abs_diff = NA_real_,
      diff_thresh = thresh,
      exceeds_thresh = NA,
      min_eig_original = NA_real_,
      min_eig_smoothed = NA_real_,
      smooth_error = paste0("eigen error: ", conditionMessage(eig0))
    ))
  }

  min_eig0 <- suppressWarnings(min(eig0, na.rm = TRUE))
  is_pd0 <- is.finite(min_eig0) && (min_eig0 > tol)

  if (is_pd0) {
    return(data.table(
      matrix = name,
      smooth_applied = FALSE,
      method = "Matrix::nearPD",
      max_abs_diff = 0,
      diff_thresh = thresh,
      exceeds_thresh = FALSE,
      min_eig_original = min_eig0,
      min_eig_smoothed = min_eig0,
      smooth_error = ""
    ))
  }

  np <- tryCatch(Matrix::nearPD(m, corr = FALSE), error = function(e) e)
  if (inherits(np, "error")) {
    return(data.table(
      matrix = name,
      smooth_applied = NA,
      method = "Matrix::nearPD",
      max_abs_diff = NA_real_,
      diff_thresh = thresh,
      exceeds_thresh = NA,
      min_eig_original = min_eig0,
      min_eig_smoothed = NA_real_,
      smooth_error = conditionMessage(np)
    ))
  }

  m2 <- as.matrix(np$mat)
  eig2 <- tryCatch(eigen(m2, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA)
  min_eig2 <- if (is.numeric(eig2)) suppressWarnings(min(eig2, na.rm = TRUE)) else NA_real_

  d <- suppressWarnings(max(abs(m - m2), na.rm = TRUE))
  exceeds <- is.finite(d) && (d > thresh)

  data.table(
    matrix = name,
    smooth_applied = TRUE,
    method = "Matrix::nearPD",
    max_abs_diff = d,
    diff_thresh = thresh,
    exceeds_thresh = exceeds,
    min_eig_original = min_eig0,
    min_eig_smoothed = min_eig2,
    smooth_error = ""
  )
}

smooth_dt <- rbind(
  smooth_report_one(mph_out$S, "S", tol = pd_tol, thresh = smooth_diff_thresh),
  smooth_report_one(mph_out$V, "V", tol = pd_tol, thresh = smooth_diff_thresh),
  fill = TRUE
)
fwrite(smooth_dt, file.path(out_dir, "smoothing_check.tsv"), sep = "\t")

smooth_error_S <- smooth_dt[matrix == "S", smooth_error]
if (length(smooth_error_S) == 1 && is.character(smooth_error_S) && nzchar(smooth_error_S)) {
  fatal_msgs <- c(fatal_msgs, paste0("S smoothing check failed: ", smooth_error_S))
}

smooth_exceeds_S <- isTRUE(smooth_dt[matrix == "S", exceeds_thresh] == TRUE)
if (smooth_exceeds_S && fail_on_smooth_exceed) {
  d <- smooth_dt[matrix == "S", max_abs_diff]
  fatal_msgs <- c(
    fatal_msgs,
    paste0(
      "S is non-PD and required smoothing is large (max |delta|=", signif(d, 4),
      " > ", smooth_diff_thresh, ")."
    )
  )
}

if (nonpd_S && fail_on_nonpd) {
  fatal_msgs <- c(
    fatal_msgs,
    paste0(
      "S is not positive definite (min eigenvalue <= ", format(pd_tol, scientific = TRUE),
      ")."
    )
  )
}

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
  if (is.null(tr)) {
    for (nm in c("trait.names", "trait_names", "traits", "var.names")) {
      if (!is.null(covstruc[[nm]])) {
        tr <- as.character(covstruc[[nm]])
        break
      }
    }
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
    "Tried: mph_out$S dimnames (rownames/colnames) and common covstruc fields; ",
    "fallback: --dict_file with column pheno_file (or trait_name).\n",
    "Your mph_out$S has no rownames/colnames and the dict_file was missing/invalid.\n",
    "Fix: pass --dict_file=/path/to/pheno_cohort_project_dict.csv"
  )
}

# A single common factor model needs at least 3 indicators to be identified under
# the default specification (F1 =~ t1 + t2 + ...). With only 2 traits, you must
# supply additional constraints via a custom MODEL_FILE.
if (length(traits) < 3 && !nzchar(model_file)) {
  msg <- paste0(
    "Only ", length(traits), " trait(s) detected.\n",
    "A 1-factor model with <3 indicators is usually under-identified under the default syntax:\n",
    "  F1 =~ trait1 + trait2\n\n",
    "Options:\n",
    "  1) Run with >=3 traits, OR\n",
    "  2) Provide MODEL_FILE with identifying constraints appropriate for a 2-indicator factor model.\n"
  )
  writeLines(msg, file.path(out_dir, "too_few_traits_FATAL.txt"))
  fatal_msgs <- c(fatal_msgs, msg)
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

# Helpful for downstream tooling; harmless if unused.
if (is.null(mph_out$trait.names)) {
  mph_out$trait.names <- traits
}

writeLines(traits, file.path(out_dir, "traits_used.txt"))

# ---- Model ----
model <- NULL
if (nzchar(model_file)) {
  if (!file.exists(model_file)) stop(paste("MODEL_FILE not found:", model_file))
  model <- paste(readLines(model_file, warn = FALSE), collapse = "\n")
} else {
  model <- paste0("F1 =~ ", paste(traits, collapse = " + "))
}

writeLines(model, file.path(out_dir, "model_used.txt"))

# ---- Fit, capturing warnings cleanly ----
warns <- character(0)
fit_error <- NULL

fit <- tryCatch({
  withCallingHandlers(
    expr = {
      GenomicSEM::usermodel(
        covstruc = mph_out,
        estimation = estimation,
        model = model,
        std.lv = std_lv
      )
    },
    warning = function(w) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(w))
      warns <<- c(warns, msg)
      invokeRestart("muffleWarning")
    }
  )
}, error = function(e) {
  fit_error <<- conditionMessage(e)
  NULL
})

warns <- unique(warns)
if (length(warns) > 0) {
  writeLines(warns, file.path(out_dir, "warnings.txt"))
} else {
  writeLines("(none)", file.path(out_dir, "warnings.txt"))
}

if (!is.null(fit_error)) {
  writeLines(fit_error, file.path(out_dir, "fit_ERROR.txt"))
  fatal_msgs <- c(fatal_msgs, paste0("usermodel() failed: ", fit_error))
}

# ---- Convergence (lavaan) ----
converged <- NA
if (!is.null(fit) && requireNamespace("lavaan", quietly = TRUE)) {
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) NA)
}

if (is.logical(converged) && length(converged) == 1 && !is.na(converged) && !converged) {
  if (fail_on_nonconverged) {
    fatal_msgs <- c(fatal_msgs, "usermodel did not converge")
  }
}

# ---- Fit stats ----
fit_stats <- NULL
if (!is.null(fit)) {
  fit_stats <- tryCatch({
    mf <- summary(fit)$modelfit
    as.data.table(t(mf), keep.rownames = "metric")
  }, error = function(e) NULL)
}

if (!is.null(fit_stats)) {
  fwrite(fit_stats, file.path(out_dir, "fit_stats.tsv"), sep = "\t")
}

# ---- Parameter estimates + Heywood flagging ----
heywood_n <- 0L
pe <- NULL
if (!is.null(fit)) {
  pe <- tryCatch({
    as.data.table(parameterEstimates(fit))
  }, error = function(e) NULL)
}

if (!is.null(pe)) {
  fwrite(pe, file.path(out_dir, "parameter_estimates.tsv"), sep = "\t")

  # Heywood heuristic: negative variances (residual or latent)
  var_rows <- pe[op == "~~" & lhs == rhs]
  heywood <- var_rows[is.finite(est) & est < 0]
  heywood_n <- nrow(heywood)
  if (heywood_n > 0) {
    fwrite(heywood, file.path(out_dir, "heywood_negative_variances.tsv"), sep = "\t")
    if (fail_on_heywood) {
      fatal_msgs <- c(fatal_msgs, paste0("Heywood case: ", heywood_n, " negative variance estimate(s)"))
    }
  }
}

if (!is.null(fit)) {
  saveRDS(fit, file.path(out_dir, "usermodel_fit.rds"))
}

# ---- Summary (single-row TSV) ----
smooth_row_S <- smooth_dt[matrix == "S"]
summary_dt <- data.table(
  n_traits = length(traits),
  nonpd_S = nonpd_S,
  smooth_applied_S = if (nrow(smooth_row_S) == 1) smooth_row_S$smooth_applied else NA,
  smooth_max_abs_diff_S = if (nrow(smooth_row_S) == 1) smooth_row_S$max_abs_diff else NA_real_,
  smooth_exceeds_thresh_S = if (nrow(smooth_row_S) == 1) smooth_row_S$exceeds_thresh else NA,
  smooth_diff_thresh = smooth_diff_thresh,
  converged = converged,
  heywood_n = heywood_n,
  fatal = (length(fatal_msgs) > 0)
)
fwrite(summary_dt, file.path(out_dir, "run_summary.tsv"), sep = "\t")

if (length(fatal_msgs) > 0) {
  fatal_msgs <- unique(fatal_msgs)
  writeLines(fatal_msgs, file.path(out_dir, "run_FATAL.txt"))
  stop(paste(fatal_msgs, collapse = "\n"))
}

message("Done. Wrote outputs to: ", out_dir)
