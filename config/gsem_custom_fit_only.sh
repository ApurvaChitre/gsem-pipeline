#!/usr/bin/env bash
# -------------------------------------------------------------------
# Custom usermodel config (2-factor example)
# Phase A only: preflight + sumstats + split + usermodel fit
# (RUN_USERGWAS=0, RUN_CREATE_DUMMY_PHENO=0)
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

# Your custom model file (2-factor). Make sure this file exists on TSCC.
# Example: config/models/model_CF_resid_both.txt
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

# Fit-only pass (no heavy array)
RUN_USERGWAS="0"
USERGWAS_ARRAY_MAX="50"

SKIP_EXISTING="1"

# -------------------------
# Report folder layout
# -------------------------
MULTIVAR_REPORT_DIR="${GSEM_DIR}/multivariate_gwas_report"
GENOME_INFO_DIR=""

# -------------------------
# Dummy pheno/data_dict (OFF here)
# -------------------------
RUN_CREATE_DUMMY_PHENO="0"
REPORT_MEASURES=""
PROCESSED_IN=""
DATA_DICT_IN=""
SOURCE_TRAIT_COL=""
SOURCE_DICT_MEASURE=""

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