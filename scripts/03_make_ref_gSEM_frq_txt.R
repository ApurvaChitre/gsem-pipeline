#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

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

frq_file <- get_arg("--frq", NA_character_)
bfile    <- get_arg("--bfile", NA_character_)  # plink prefix, optional
bim_file <- get_arg("--bim", NA_character_)
out_file <- get_arg("--out", NA_character_)

if (is.na(frq_file) || frq_file == "") {
  stop("Missing required arg: --frq=/path/to/ref_gSEM.frq")
}
if (is.na(out_file) || out_file == "") {
  stop("Missing required arg: --out=/path/to/ref_gSEM_frq.txt")
}

# If user gave --bfile but not --bim, infer bim path
if ((is.na(bim_file) || bim_file == "") && !is.na(bfile) && bfile != "") {
  bim_file <- paste0(bfile, ".bim")
}

if (!file.exists(frq_file)) stop("frq file not found: ", frq_file)

frq <- fread(frq_file)

# Expect PLINK .frq columns
need_frq <- c("CHR", "SNP", "A1", "A2", "MAF")
miss_frq <- setdiff(need_frq, names(frq))
if (length(miss_frq) > 0) {
  stop(".frq missing columns: ", paste(miss_frq, collapse = ", "))
}

# Get BP either from BIM (preferred) or from parsing SNP if it is CHR:BP
bp_dt <- NULL
if (!is.na(bim_file) && bim_file != "") {
  if (!file.exists(bim_file)) stop("bim file not found: ", bim_file)
  bim <- fread(bim_file, header = FALSE)
  if (ncol(bim) < 4) stop("Unexpected .bim format (need at least 4 columns): ", bim_file)
  setnames(bim, c("CHR_bim", "SNP", "CM", "BP", "A1_bim", "A2_bim")[1:ncol(bim)])
  bp_dt <- bim[, .(SNP, BP = as.integer(BP))]
}

if (!is.null(bp_dt)) {
  dt <- merge(frq[, .(SNP, CHR, MAF, A1, A2)], bp_dt, by = "SNP", all.x = TRUE)
} else {
  # Fallback: parse BP from SNP like "1:735254".
  parts <- tstrsplit(frq$SNP, ":", fixed = TRUE)
  if (length(parts) < 2) {
    stop("No --bim provided and SNP IDs don't look like CHR:BP (e.g. 1:735254). Provide --bim or --bfile.")
  }
  dt <- frq[, .(
    SNP,
    CHR = as.integer(CHR),
    BP  = suppressWarnings(as.integer(parts[[2]])),
    MAF,
    A1,
    A2
  )]
}

if (any(is.na(dt$BP))) {
  bad <- dt[is.na(BP), head(SNP, 10)]
  stop("BP is NA for some SNPs (showing up to 10): ", paste(bad, collapse = ", "),
       "\nProvide a .bim via --bim=... or --bfile=... to resolve positions.")
}

# Add chr prefix to SNP IDs if not already present
# Example: "1:735254" -> "chr1:735254"; "rs123" -> "chrrs123" is NOT desired,
# so only prefix when SNP looks like coordinate (contains ':') or starts with digit.
looks_coord <- grepl(":", dt$SNP) | grepl("^[0-9]", dt$SNP)
SNP_out <- dt$SNP
SNP_out[looks_coord & !grepl("^chr", SNP_out)] <- paste0("chr", SNP_out[looks_coord & !grepl("^chr", SNP_out)])

out <- dt[, .(
  SNP = SNP_out,
  CHR = as.integer(CHR),
  BP  = as.integer(BP),
  MAF = as.numeric(MAF),
  A1,
  A2
)]

# Write as tab-delimited
fwrite(out, out_file, sep = "\t", quote = FALSE)
cat("Wrote:", out_file, "\n")


