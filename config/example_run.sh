#!/usr/bin/env bash
# -------------------------------------------------------------------
# GenomicSEM usermodel (1-factor) pipeline config (TSCC / Slurm)
#
# This file is a concrete "test run" config using Apurva's paper_2026 paths.
# Coworkers should copy it and change MPH_DIR / GSEM_DIR / UNIV_GWAS_DIR / GENO_BFILE
# to their own locations.
#
# Run:
#   bash bin/submit_debug_env.sh --config config/example_run.sh
#   bash bin/submit_gsem_pipeline.sh --config config/example_run.sh
# -------------------------------------------------------------------

# -------------------------
# REQUIRED INPUTS
# -------------------------

# MPH pipeline output directory (must contain gsem/ + pheno/)
MPH_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/mph_TEST_run1"

# Optional: explicit trait dictionary file (pheno_cohort_project_dict.csv).
# If empty, the pipeline will try common MPH layouts under MPH_DIR.
DICT_FILE=""

# Where to write GenomicSEM outputs for this run
GSEM_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/gsem_TEST_run"

# Directory containing univariate GWAS *.mlma files (searched recursively)
UNIV_GWAS_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/individual_projects/results/gwas"

# Genotype prefix for plink (must have .bed/.bim/.fam)
GENO_BFILE="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/individual_projects/genotypes/genotypes"

# Autosomes
CHR_RANGE="1-20"

# Pattern used to select formatted sumstats files for GenomicSEM::sumstats()
SUMSTATS_PATTERN="_formatted.txt$"


# -------------------------
# OPTIONAL: conda environment
# -------------------------

CONDA_ENV="genomic_sem"

# Optional: path to conda.sh for non-interactive shells.
# If conda activation fails inside Slurm, set this explicitly:
#   CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
CONDA_SH=""

# If you want to force using a specific plink binary, set it here
PLINK_BIN="plink"

# -------------------------
# MODEL
# -------------------------

# Used for naming output folders under results/multivariate_gwas/
MODEL_NAME="usermodel_1factor"

# Optional: text file containing lavaan model syntax.
# If empty, the pipeline auto-builds:
#   F1 =~ trait1 + trait2 + ... + traitN
MODEL_FILE=""

# -------------------------
# SPLITTING / PIPELINE FLOW
# -------------------------

# -------------------------
# STEP 0: PREFLIGHT (recommended)
# -------------------------

# If 1, run a quick sanity check on the MPH/GenomicSEM input matrices (S/V)
# *before* we spend time building sumstats + splitting SNP chunks.
RUN_PREFLIGHT="1"

# Policy (matches common GenomicSEM practice):
#   - If S is non-PD but the implied smoothing is small, proceed (warn).
#   - If smoothing is large (default threshold 0.025), stop early.
#   - Stop early for obviously pathological model fits (non-convergence, Heywood cases) so you
#     don't waste compute on the userGWAS array.

# Hard-stop immediately on non-PD S based only on eigenvalues (usually NOT needed).
# Leave this 0 unless you want the strictest behavior.
PREFLIGHT_FAIL_ON_NONPD="0"

# If 1 (default), stop when the max absolute change required by nearPD smoothing exceeds threshold.
PREFLIGHT_FAIL_ON_SMOOTH_EXCEED="1"

# Smoothing threshold used by GenomicSEM warnings in many workflows.
PREFLIGHT_SMOOTH_DIFF_THRESH="0.025"

# Stop on model pathologies
PREFLIGHT_FAIL_ON_HEYWOOD="1"
PREFLIGHT_FAIL_ON_NONCONVERGED="1"

# Tolerance for PD check (min eigenvalue must be > this threshold to count as PD)
PREFLIGHT_PD_TOL="-1e-8"

# SNPs per chunk for userGWAS
# (Bigger = fewer array tasks; smaller = less memory per task)
SPLIT_SIZE="100000"

# Optional: limit how many userGWAS array tasks run concurrently.
# Set empty or 0 for no throttle.
USERGWAS_ARRAY_MAX="50"


# If 1, run usermodel fit once before the userGWAS array
RUN_USERMODEL_FIT="1"

# IMPORTANT: userGWAS is the heavy step.
# Recommended workflow:
#   1) RUN_USERGWAS="0"  -> prep + split + model-fit ONLY (stops before userGWAS)
#   2) inspect ${GSEM_DIR}/results/model_fit/
#   3) set RUN_USERGWAS="1" -> re-run submit_gsem_pipeline.sh to launch userGWAS
RUN_USERGWAS="0"

# If 1, skip steps that already have a .done sentinel file
SKIP_EXISTING="1"

# -------------------------
# Downstream reporting folder layout
# -------------------------

# This directory will be created and populated with compiled *.mlma files.
# It matches the folder layout expected by multivariate_gwas_report.ipynb.
MULTIVAR_REPORT_DIR="${GSEM_DIR}/multivariate_gwas_report"

# Optional: if you have a genome_info/ directory from the univariate pipeline,
# provide it and we will copy it into MULTIVAR_REPORT_DIR/genome_info/
GENOME_INFO_DIR=""

# Optional: genotype dir to symlink into MULTIVAR_REPORT_DIR/genotypes/
REPORT_GENOTYPES_DIR="/tscc/projects/ps-palmer/apurva/locomotor/paper_2026/individual_projects/genotypes"

# -------------------------
# Optional: create dummy processed_data_ready + data dict for report pipeline
# -------------------------
RUN_CREATE_DUMMY_PHENO="0"
REPORT_MEASURES=""
PROCESSED_IN=""
DATA_DICT_IN=""
SOURCE_TRAIT_COL=""
SOURCE_DICT_MEASURE=""

# -------------------------
# Slurm defaults (Apurva's working TSCC condo allocation)
# -------------------------
SLURM_ACCOUNT="csd795"
SLURM_PARTITION="condo"
SLURM_QOS="condo"

# -------------------------
# Slurm resources
# These mirror the "manual" srun/sbatch settings you were using in paper_2026.
# -------------------------

# 02_make_gwas_sumstats_final.R (merge/format univariate GWAS)
MAKE_SUMSTATS_TIME="1:00:00"; MAKE_SUMSTATS_CPUS="1"; MAKE_SUMSTATS_MEM="32G"

# 03_make_ref_gsem (plink freq + convert)
MAKE_REF_TIME="1:00:00";     MAKE_REF_CPUS="1";     MAKE_REF_MEM="16G"

# 04_sumstats_func_gsem.R (GenomicSEM::sumstats)
SUMSTATS_TIME="1:00:00";     SUMSTATS_CPUS="1";     SUMSTATS_MEM="16G"

# split sumstats into chunks
SPLIT_TIME="3:00:00";        SPLIT_CPUS="8";        SPLIT_MEM="60G"

# preflight (cheap)
PREFLIGHT_TIME="0:30:00";    PREFLIGHT_CPUS="1";    PREFLIGHT_MEM="8G"

# usermodel fit (cheap gate)
USERMODEL_TIME="0:30:00";    USERMODEL_CPUS="1";    USERMODEL_MEM="32G"

# userGWAS array (heavy)
USERGWAS_TIME="12:00:00";    USERGWAS_CPUS="16";    USERGWAS_MEM="64G"

# compile chunk outputs into per-chr MLMA
COMPILE_TIME="1:00:00";      COMPILE_CPUS="1";      COMPILE_MEM="32G"

# assemble multivariate_gwas_report/ folder
REPORT_TIME="0:30:00";       REPORT_CPUS="1";       REPORT_MEM="8G"
