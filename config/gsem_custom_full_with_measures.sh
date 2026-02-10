#!/usr/bin/env bash
# -------------------------------------------------------------------
# Custom usermodel config (2-factor example)
# Full pipeline: includes userGWAS + compile + report + dummy pheno/data_dict
# -------------------------------------------------------------------

# -------------------------
# REQUIRED INPUTS
# -------------------------
MPH_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph_TEST_run1"
DICT_FILE=""

GSEM_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/gsem_TEST_run"
UNIV_GWAS_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/individual_projects/results/gwas"
GENO_BFILE="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/individual_projects/genotypes/genotypes"

CHR_RANGE="1-20"
SUMSTATS_PATTERN="_formatted.txt$"

# -------------------------
# Conda environment
# -------------------------
CONDA_ENV="genomic_sem"
CONDA_SH=""
PLINK_BIN="plink"

# -------------------------
# MODEL
# -------------------------
MODEL_NAME="cf_resid_both"
MODEL_FILE="config/models/model_CF_resid_both.txt"

# -------------------------
# PREFLIGHT / FLOW
# -------------------------
SPLIT_SIZE="100000"

RUN_PREFLIGHT="1"
PREFLIGHT_FAIL_ON_NONPD="0"
PREFLIGHT_FAIL_ON_SMOOTH_EXCEED="1"
PREFLIGHT_SMOOTH_DIFF_THRESH="0.025"
PREFLIGHT_FAIL_ON_HEYWOOD="1"
PREFLIGHT_FAIL_ON_NONCONVERGED="1"
PREFLIGHT_PD_TOL="-1e-8"

RUN_USERMODEL_FIT="1"

# Full run (heavy)
RUN_USERGWAS="1"
USERGWAS_ARRAY_MAX="50"

SKIP_EXISTING="1"

# -------------------------
# Report folder layout
# -------------------------
MULTIVAR_REPORT_DIR="${GSEM_DIR}/multivariate_gwas_report"
GENOME_INFO_DIR=""

# -------------------------
# Dummy pheno/data_dict (ON)
# IMPORTANT:
# Measures are WITHOUT regressedlr_ prefix.
# Must match MLMA basenames (lowercase).
# For 2-factor model (F1/F2):
#   cf_resid_both_usergwas_f1
#   cf_resid_both_usergwas_f2
# -------------------------
RUN_CREATE_DUMMY_PHENO="1"
#REPORT_MEASURES="cf_resid_both_usergwas_f1,cf_resid_both_usergwas_f2"
REPORT_MEASURES="cf_resid_both_usergwas_f1"

# These are your "mega" inputs to copy from (edit as needed)
PROCESSED_IN="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mega/processed_data_ready.csv"
DATA_DICT_IN="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mega/data_dict_mega.csv"
SOURCE_TRAIT_COL="regressedlr_locomotor_mega"
SOURCE_DICT_MEASURE="locomotor_mega"

# -------------------------
# Slurm defaults
# -------------------------
SLURM_ACCOUNT="csd795"
SLURM_PARTITION="condo"
SLURM_QOS="condo"

# -------------------------
# Resources
# -------------------------
PREFLIGHT_TIME="0:30:00";      PREFLIGHT_CPUS="1";      PREFLIGHT_MEM="8G"
MAKE_SUMSTATS_TIME="1:00:00";  MAKE_SUMSTATS_CPUS="1";  MAKE_SUMSTATS_MEM="32G"
MAKE_REF_TIME="1:00:00";       MAKE_REF_CPUS="1";       MAKE_REF_MEM="16G"
SUMSTATS_TIME="1:00:00";       SUMSTATS_CPUS="1";       SUMSTATS_MEM="16G"
SPLIT_TIME="3:00:00";          SPLIT_CPUS="8";          SPLIT_MEM="60G"
USERMODEL_TIME="0:30:00";      USERMODEL_CPUS="1";      USERMODEL_MEM="32G"
USERGWAS_TIME="12:00:00";      USERGWAS_CPUS="16";      USERGWAS_MEM="64G"
COMPILE_TIME="1:00:00";        COMPILE_CPUS="1";        COMPILE_MEM="32G"
REPORT_TIME="0:30:00";         REPORT_CPUS="1";         REPORT_MEM="8G"