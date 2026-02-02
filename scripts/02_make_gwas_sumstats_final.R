#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
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

has_flag <- function(flag) {
  any(args == flag) || any(startsWith(args, paste0(flag, '=')))
}

if (has_flag('--help') || has_flag('-h')) {
  cat(paste0(
    'Make final GWAS summary stats (merge chr files + format)\n\n',
    'Usage:\n',
    '  Rscript make_gwas_sumstats_final.R ',
    '--gwas_dir=/path/to/gwas ',
    '--dict_file=/path/to/pheno_cohort_project_dict.csv ',
    '--out_dir=/path/to/sum_stats_final ',
    '[--chr_min=1 --chr_max=20 --snp_prefix=chr]\n\n',
    'Notes:\n',
    '  - Expects per-chrom files named like: regressedlr_<trait>_chrgwas<CHR>.mlma\n',
    '  - Uses N from the dict_file for each trait and writes one formatted file per trait.\n'
  ))
  quit(status = 0)
}

# Defaults (match your existing pipeline)
gwas_dir  <- get_arg('--gwas_dir',  '/tscc/projects/ps-palmer/apurva/locomotor/individual_projects_seq/results/gwas')
dict_file <- get_arg('--dict_file', '/tscc/projects/ps-palmer/apurva/locomotor/mph_ccc_kalivas/individual_projects_seq/pheno_cohort_project_dict.csv')
out_dir   <- get_arg('--out_dir',   '/tscc/projects/ps-palmer/apurva/locomotor/mph_ccc_kalivas/individual_projects_seq/sum_stats_final')

chr_min <- as.integer(get_arg('--chr_min', '1'))
chr_max <- as.integer(get_arg('--chr_max', '20'))
chr_vec <- chr_min:chr_max

snp_prefix <- get_arg('--snp_prefix', 'chr')  # set to empty to disable prefix

# Validate inputs
if (is.na(gwas_dir) || !dir.exists(gwas_dir)) stop('gwas_dir does not exist: ', gwas_dir)
if (is.na(dict_file) || !file.exists(dict_file)) stop('dict_file does not exist: ', dict_file)
if (is.na(out_dir) || nchar(out_dir) == 0) stop('out_dir is empty')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trait_dict <- fread(dict_file)

# Required columns
if (!('pheno_file' %in% names(trait_dict))) stop('dict_file must contain a pheno_file column: ', dict_file)
if (!('N' %in% names(trait_dict))) stop('dict_file must contain an N column (sample size): ', dict_file)

# Clean trait names: remove optional regressedlr_ prefix + .txt suffix
trait_dict[, trait_name := pheno_file]
trait_dict[, trait_name := str_trim(trait_name)]
trait_dict[, trait_name := gsub('\\.txt$', '', trait_name)]
trait_dict[, trait_name := gsub('^regressedlr_', '', trait_name)]

# Make sure we iterate only over unique traits (defensive)
trait_dict <- unique(trait_dict, by = 'trait_name')

required_mlma_cols <- c('Chr', 'SNP', 'A1', 'A2', 'b', 'se', 'p')

n_done <- 0L
n_missing <- 0L

for (i in 1:nrow(trait_dict)) {
  trait_name  <- trait_dict$trait_name[i]
  sample_size <- trait_dict$N[i]

  if (is.na(trait_name) || trait_name == '') {
    message('Skipping blank trait_name at row ', i)
    next
  }

  # Build expected file paths in chr order (avoids regex edge cases)
  candidate_files <- file.path(
    gwas_dir,
    paste0('regressedlr_', trait_name, '_chrgwas', chr_vec, '.mlma')
  )
  chrom_files <- candidate_files[file.exists(candidate_files)]

  if (length(chrom_files) == 0) {
    message('No chromosome files found for: ', trait_name)
    n_missing <- n_missing + 1L
    next
  }

  # Read + concatenate
  gwas_data <- rbindlist(lapply(chrom_files, fread), use.names = TRUE, fill = TRUE)

  # Basic schema check
  missing_cols <- setdiff(required_mlma_cols, names(gwas_data))
  if (length(missing_cols) > 0) {
    message('Skipping ', trait_name, ' (missing columns in mlma): ', paste(missing_cols, collapse = ', '))
    n_missing <- n_missing + 1L
    next
  }

  # Keep only requested chromosomes  (Chr can be character in some outputs)
  gwas_data[, Chr_int := suppressWarnings(as.integer(Chr))]
  gwas_data <- gwas_data[Chr_int %in% chr_vec]

  # Format final sumstats
  formatted <- gwas_data[, .(
    SNP = if (nchar(snp_prefix) > 0) paste0(snp_prefix, SNP) else as.character(SNP),
    A1 = toupper(A1),
    A2 = toupper(A2),
    effect = b,
    SE = se,
    P = p,
    N = sample_size
  )]

  # Drop missing
  formatted <- formatted[!is.na(SNP) & !is.na(effect) & !is.na(SE) & !is.na(P)]

  # Write
  output_file <- file.path(out_dir, paste0(trait_name, '_formatted.txt'))
  fwrite(formatted, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)

  n_done <- n_done + 1L
  if (length(chrom_files) < length(chr_vec)) {
    # Light warning if some chr files are absent (still wrote output)
    message('Wrote ', basename(output_file), ' (', length(chrom_files), '/', length(chr_vec), ' chr files present)')
  } else {
    message('Wrote ', basename(output_file), ' (', length(chrom_files), ' chr files)')
  }
}

message('--- Summary ---')
message('Traits written: ', n_done)
message('Traits skipped/missing: ', n_missing)
message('Output directory: ', out_dir)


