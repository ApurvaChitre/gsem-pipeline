#!/usr/bin/env bash
#SBATCH --job-name=gsem_make_sumstats
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

# ----- Conda activation (recommended: provides R packages like data.table) -----
if [[ -n "${CONDA_ENV:-}" ]]; then
  if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
  else
    # On TSCC, conda is usually initialized in ~/.bashrc
    [[ -f "${HOME}/.bashrc" ]] && source "${HOME}/.bashrc"
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
  fi
  conda activate "${CONDA_ENV}"
  echo "[env] Activated conda env: ${CONDA_ENV}"
  echo "[env] Rscript: $(command -v Rscript || true)"
  echo "[env] plink:   $(command -v plink || true)"
fi

mkdir -p "${GSEM_DIR}/sum_stats_final" "${GSEM_DIR}/logs"

DONE_FILE="${GSEM_DIR}/sum_stats_final/.make_sumstats.done"
if [[ "${SKIP_EXISTING}" == "1" && -f "${DONE_FILE}" ]]; then
  echo "[make_sumstats] SKIP_EXISTING=1 and found ${DONE_FILE}; skipping."
  exit 0
fi

# Dictionary (trait list + N). Prefer explicit DICT_FILE in config; fallback to common MPH layouts.
DICT_F="${DICT_FILE:-}"
if [[ -z "${DICT_F}" ]]; then
  DICT_F="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"
fi
if [[ ! -f "${DICT_F}" ]]; then
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

if [[ ! -f "${DICT_F}" ]]; then
  echo "ERROR: pheno_cohort_project_dict.csv not found." >&2
  echo "  Looked for: \"${DICT_F}\" and common MPH defaults under MPH_DIR=${MPH_DIR}" >&2
  echo "  Fix: set DICT_FILE=\"/path/to/pheno_cohort_project_dict.csv\" in your config." >&2
  exit 2
fi

echo "[make_sumstats] UNIV_GWAS_DIR=${UNIV_GWAS_DIR}"
echo "[make_sumstats] DICT_FILE=${DICT_F}"
echo "[make_sumstats] out_dir=${GSEM_DIR}/sum_stats_final"

Rscript "${REPO_DIR}/scripts/02_make_gwas_sumstats_final.R" \
  --gwas_dir="${UNIV_GWAS_DIR}" \
  --dict_file="${DICT_F}" \
  --out_dir="${GSEM_DIR}/sum_stats_final"

touch "${DONE_FILE}"
echo "[make_sumstats] DONE"
