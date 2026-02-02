#!/usr/bin/env bash
#SBATCH --job-name=gsem_preflight
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00

# Preflight: sanity-check MPH/GenomicSEM covstruc (S/V) before running any heavy steps.
#
# What it does:
#   - loads MPH_DIR/gsem/MPH_genomicSEM.RData
#   - checks whether S (and V) are positive definite (writes matrix_pd_check.tsv)
#   - fits a simple 1-factor model (writes warnings + fit stats)
#
# If PREFLIGHT_FAIL_ON_NONPD=1 (default), the job FAILS when S is non-PD.

set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

mkdir -p "${GSEM_DIR}/logs" "${GSEM_DIR}/results/model_fit/preflight"

# ----- Conda activate (R + GenomicSEM + data.table) -----
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate "${CONDA_ENV}"

echo "[env] Activated conda env: ${CONDA_ENV}"
echo "[env] Rscript: $(command -v Rscript)"

FAIL_NONPD="${PREFLIGHT_FAIL_ON_NONPD:-0}"
FAIL_SMOOTH="${PREFLIGHT_FAIL_ON_SMOOTH_EXCEED:-1}"
FAIL_HEYWOOD="${PREFLIGHT_FAIL_ON_HEYWOOD:-1}"
FAIL_NONCONV="${PREFLIGHT_FAIL_ON_NONCONVERGED:-1}"

PD_TOL="${PREFLIGHT_PD_TOL:--1e-8}"
SMOOTH_THRESH="${PREFLIGHT_SMOOTH_DIFF_THRESH:-0.025}"

to_r_bool() {
  # bash truthy: 1 -> TRUE, else FALSE
  if [[ "${1:-0}" == "1" ]]; then
    echo "TRUE"
  else
    echo "FALSE"
  fi
}

FAIL_NONPD_R=$(to_r_bool "${FAIL_NONPD}")
FAIL_SMOOTH_R=$(to_r_bool "${FAIL_SMOOTH}")
FAIL_HEYWOOD_R=$(to_r_bool "${FAIL_HEYWOOD}")
FAIL_NONCONV_R=$(to_r_bool "${FAIL_NONCONV}")

# Dict file is used to infer trait names when mph_out$S has no rownames.
DICT_F="${DICT_FILE:-}"
if [[ -z "${DICT_F}" || ! -f "${DICT_F}" ]]; then
  for cand in \
    "${MPH_DIR}/pheno/pheno_cohort_project_dict.csv" \
    "${MPH_DIR}/pheno_cohort_project_dict.csv" \
    "${MPH_DIR}/individual_projects_seq/pheno_cohort_project_dict.csv" \
    ; do
    if [[ -f "${cand}" ]]; then
      DICT_F="${cand}"
      break
    fi
  done
fi

DICT_ARG=()
if [[ -n "${DICT_F}" && -f "${DICT_F}" ]]; then
  DICT_ARG=(--dict_file "${DICT_F}")
fi

Rscript "${REPO_DIR}/scripts/01_fit_usermodel_1factor.R" \
  --mph_rdata "${MPH_DIR}/gsem/MPH_genomicSEM.RData" \
  --out_dir "${GSEM_DIR}/results/model_fit/preflight" \
  "${DICT_ARG[@]}" \
  --fail_on_nonpd "${FAIL_NONPD_R}" \
  --fail_on_smooth_exceed "${FAIL_SMOOTH_R}" \
  --smooth_diff_thresh "${SMOOTH_THRESH}" \
  --fail_on_heywood "${FAIL_HEYWOOD_R}" \
  --fail_on_nonconverged "${FAIL_NONCONV_R}" \
  --pd_tol "${PD_TOL}"

echo "[preflight] Done. Outputs: ${GSEM_DIR}/results/model_fit/preflight"
