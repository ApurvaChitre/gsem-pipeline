#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------
# 0) Args: --input_dir=... --out_dir=... [--prefix=...]
# ----------------------------------------------------------------------------

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

input_dir <- get_arg("--input_dir", NA_character_)
out_dir   <- get_arg("--out_dir",   NA_character_)
prefix    <- get_arg("--prefix",    "regressedlr_gsem_results")

if (is.na(input_dir) || input_dir == "") stop("Missing required --input_dir=/path")
if (is.na(out_dir)   || out_dir   == "") stop("Missing required --out_dir=/path")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_slug <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

suppressPackageStartupMessages(library(dplyr))

cat("Input dir:", input_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Prefix:", prefix, "\n")

# ----------------------------------------------------------------------------
# 1) List all .RData/.rdata in input_dir and group by model (chunk-aware)
#    If filename matches: <model>_chunk_<N>.RData => grouped and sorted by chunk
#    Else: treated as its own model group (one file)
# ----------------------------------------------------------------------------

rdata_files <- list.files(input_dir, pattern = "\\.[Rr]Data$", full.names = TRUE)
if (length(rdata_files) == 0) stop("No .RData files found in: ", input_dir)

rx <- "^(.*)_chunk_([0-9]+)\\.[Rr]Data$"
m <- regexec(rx, basename(rdata_files))
parts <- regmatches(basename(rdata_files), m)
is_chunked <- lengths(parts) == 3

# Map files -> model_pat + chunk
file_map <- data.frame(
  file      = rdata_files,
  model_pat = ifelse(is_chunked,
                     vapply(parts[is_chunked], `[[`, "", 2),
                     tools::file_path_sans_ext(basename(rdata_files[!is_chunked]))),
  chunk     = NA_integer_,
  stringsAsFactors = FALSE
)

# Fill chunk numbers where chunked
file_map$chunk[is_chunked] <- as.integer(vapply(parts[is_chunked], `[[`, "", 3))

# For non-chunked files, model_pat was filled only for those rows via ifelse above,
# but because rdata_files[!is_chunked] is shorter, fix model_pat for all rows explicitly:
file_map$model_pat[!is_chunked] <- tools::file_path_sans_ext(basename(file_map$file[!is_chunked]))

model_patterns <- sort(unique(file_map$model_pat))
cat("Discovered model groups:", paste(model_patterns, collapse = ", "), "\n")

# ----------------------------------------------------------------------------
# 2) Your parse functions (kept, with minimal robustness improvements)
# ----------------------------------------------------------------------------

# A) parse_commonfactor(): for commonfactor-like data frames/lists
parse_commonfactor <- function(cf_obj) {
  if (is.matrix(cf_obj)) cf_obj <- as.data.frame(cf_obj)
  
  if (is.list(cf_obj) && !is.data.frame(cf_obj)) {
    cf_obj <- do.call(rbind, cf_obj)
  }
  
  # tolerate SE/se instead of se_c
  if (!"se_c" %in% names(cf_obj)) {
    if ("SE" %in% names(cf_obj)) cf_obj$se_c <- cf_obj$SE
    if ("se" %in% names(cf_obj)) cf_obj$se_c <- cf_obj$se
  }
  
  keep_cols <- c("CHR", "BP", "SNP", "MAF", "A1", "A2",
                 "est", "se_c", "Z_Estimate", "Pval_Estimate")
  miss <- setdiff(keep_cols, names(cf_obj))
  if (length(miss)) stop("parse_commonfactor(): missing cols: ", paste(miss, collapse = ", "))
  
  df <- cf_obj[, keep_cols]
  names(df)[names(df) == "se_c"]          <- "SE"
  names(df)[names(df) == "Pval_Estimate"] <- "pval"
  df$factor <- "F1_common"
  df
}

# B) parse_user_model(): for factor_out-like data frames/lists
parse_user_model <- function(factor_list) {
  # factor_list is typically a list of small data frames
  if (is.matrix(factor_list)) factor_list <- as.data.frame(factor_list)
  if (is.data.frame(factor_list)) factor_list <- list(factor_list)
  
  df_list <- lapply(factor_list, function(df_chunk) {
    if (is.matrix(df_chunk)) df_chunk <- as.data.frame(df_chunk)
    if ("se" %in% names(df_chunk) && !"SE" %in% names(df_chunk)) {
      names(df_chunk)[names(df_chunk) == "se"] <- "SE"
    }
    keep <- (df_chunk$op == "~" & df_chunk$rhs == "SNP")
    df_chunk[keep, ]
  })
  
  df <- do.call(rbind, df_list)
  keep_cols <- c("CHR", "BP", "SNP", "MAF", "A1", "A2",
                 "lhs", "est", "SE", "Z_Estimate", "Pval_Estimate")
  miss <- setdiff(keep_cols, names(df))
  if (length(miss)) stop("parse_user_model(): missing cols: ", paste(miss, collapse = ", "))
  
  df <- df[, keep_cols]
  names(df)[names(df) == "lhs"]           <- "factor"
  names(df)[names(df) == "Pval_Estimate"] <- "pval"
  df
}

# ----------------------------------------------------------------------------
# 3) Structure-based detection (no dependence on object names)
# ----------------------------------------------------------------------------

looks_like_commonfactor <- function(x) {
  req <- c("CHR","BP","SNP","MAF","A1","A2","est","Pval_Estimate")
  se_ok <- function(nm) any(c("se_c","SE","se") %in% nm)
  
  if (is.matrix(x)) x <- as.data.frame(x)
  
  if (is.data.frame(x)) {
    nm <- names(x)
    return(all(req %in% nm) && se_ok(nm))
  }
  
  if (is.list(x) && length(x) > 0 &&
      all(vapply(x, function(z) is.data.frame(z) || is.matrix(z), TRUE))) {
    z <- x[[1]]
    if (is.matrix(z)) z <- as.data.frame(z)
    nm <- names(z)
    return(all(req %in% nm) && se_ok(nm))
  }
  
  FALSE
}

looks_like_factor_out <- function(x) {
  req <- c("op","rhs","lhs","CHR","BP","SNP","MAF","A1","A2","est","Pval_Estimate")
  se_ok <- function(nm) any(c("SE","se") %in% nm)
  
  if (is.matrix(x)) x <- as.data.frame(x)
  
  if (is.data.frame(x)) {
    nm <- names(x)
    return(all(req %in% nm) && se_ok(nm))
  }
  
  if (is.list(x) && length(x) > 0 &&
      all(vapply(x, function(z) is.data.frame(z) || is.matrix(z), TRUE))) {
    z <- x[[1]]
    if (is.matrix(z)) z <- as.data.frame(z)
    nm <- names(z)
    return(all(req %in% nm) && se_ok(nm))
  }
  
  FALSE
}

try_parse_any_loaded_object <- function(obj) {
  # Prefer factor_out if it matches (more distinctive)
  if (looks_like_factor_out(obj)) {
    return(tryCatch(parse_user_model(obj), error = function(e) NULL))
  }
  if (looks_like_commonfactor(obj)) {
    return(tryCatch(parse_commonfactor(obj), error = function(e) NULL))
  }
  NULL
}

# ----------------------------------------------------------------------------
# 4) For each model group: load all files, parse, combine, write MLMA per factor×chr
# ----------------------------------------------------------------------------

for (model_pat in model_patterns) {
  cat("\n----------------------------------------\n")
  cat("Processing model group:", model_pat, "\n")
  
  this_map <- file_map[file_map$model_pat == model_pat, , drop = FALSE]
  
  # Sort: chunked by chunk; non-chunked keep filename order
  if (all(is.na(this_map$chunk))) {
    file_vec <- this_map$file
  } else {
    # chunked files first in ascending chunk; any NA chunks last
    ord <- order(is.na(this_map$chunk), this_map$chunk, basename(this_map$file))
    file_vec <- this_map$file[ord]
    
    # warn if missing chunks (only meaningful when chunked)
    chunk_vals <- sort(unique(this_map$chunk[!is.na(this_map$chunk)]))
    if (length(chunk_vals) > 0) {
      expected <- seq_len(max(chunk_vals))
      missing <- setdiff(expected, chunk_vals)
      if (length(missing) > 0) {
        cat("WARNING: missing chunk numbers for", model_pat, ":\n")
        print(missing)
      }
    }
  }
  
  all_chunks <- list()
  
  for (f in file_vec) {
    cat("  Loading:", basename(f), "...\n")
    
    e <- new.env(parent = emptyenv())
    objs_loaded <- load(f, envir = e)
    
    df_parsed <- NULL
    for (obj_name in objs_loaded) {
      obj <- get(obj_name, envir = e)
      df_parsed <- try_parse_any_loaded_object(obj)
      if (!is.null(df_parsed)) {
        cat("    Parsed object by structure:", obj_name, "\n")
        break
      }
    }
    
    if (is.null(df_parsed)) {
      cat("    ERROR: No recognized object found in", basename(f), "- skipping.\n")
      next
    }
    
    all_chunks[[length(all_chunks) + 1]] <- df_parsed
    rm(e)
  }
  
  if (length(all_chunks) == 0) {
    cat("  No valid chunks parsed for model group", model_pat, "- skipping.\n")
    next
  }
  
  model_res <- do.call(rbind, all_chunks)
  
  # Validate minimal columns needed for MLMA conversion
  needed <- c("CHR","BP","SNP","A1","A2","MAF","est","SE","pval","factor")
  miss <- setdiff(needed, names(model_res))
  if (length(miss)) stop("After parsing, missing cols: ", paste(miss, collapse = ", "))
  
  # Convert to MLMA schema in memory; keep factor only for splitting
  df <- model_res %>%
    rename(
      Chr   = CHR,
      bp    = BP,
      Freq  = MAF,
      b     = est,
      se    = SE,
      p     = pval,
      trait = factor
    ) %>%
    mutate(SNP = sub("^chr", "", SNP)) %>%  # strip chr prefix if present
    select(Chr, bp, SNP, A1, A2, Freq, b, se, p, trait)
  
  safe_model <- safe_slug(model_pat)
  
  # Write per factor x chr MLMA (trait NOT inside file)
  factors <- unique(df$trait)
  cat("  Factors detected:", paste(factors, collapse = ", "), "\n")
  
  for (fac in factors) {
    df_fac <- df %>% filter(trait == fac)
    safe_fac <- safe_slug(fac)
    
    chrs <- sort(unique(df_fac$Chr))
    for (chr in chrs) {
      df_chr <- df_fac %>%
        filter(Chr == chr) %>%
        select(Chr, bp, SNP, A1, A2, Freq, b, se, p)
      
      out_file <- file.path(
        out_dir,
        sprintf("%s_%s_%s_chrgwas%d.mlma", prefix, safe_model, safe_fac, chr)
      )
      
      write.table(df_chr, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  cat("  Done writing MLMA files for model group:", model_pat, "\n")
}

cat("\nAll done!\n")

