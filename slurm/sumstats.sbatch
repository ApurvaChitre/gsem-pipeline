#!/usr/bin/env bash
#SBATCH --job-name=gsem_sumstats
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

mkdir -p "${GSEM_DIR}/sum_stats_final" "${GSEM_DIR}/logs"

OUT_RDATA="${GSEM_DIR}/sum_stats_final/mySumstatsGSEM.RData"
DONE_FILE="${GSEM_DIR}/sum_stats_final/.sumstats.done"

if [[ "${SKIP_EXISTING}" == "1" && -f "${DONE_FILE}" && -s "${OUT_RDATA}" ]]; then
  echo "[sumstats] DONE_FILE exists and mySumstatsGSEM.RData non-empty; skipping."
  exit 0
fi

# ----- Conda activate -----
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate "${CONDA_ENV}"

# Default: only use formatted sumstats files
PATTERN="${SUMSTATS_PATTERN:-_formatted.txt$}"

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
  echo "  Looked for common MPH defaults under MPH_DIR=${MPH_DIR}" >&2
  echo "  Fix: set DICT_FILE=\"/path/to/pheno_cohort_project_dict.csv\" in your config." >&2
  exit 2
fi


Rscript "${REPO_DIR}/scripts/04_sumstats_func_gsem.R" \
  --base_dir="${GSEM_DIR}" \
  --dict_file="${DICT_F}" \
  --sumstats_dir="${GSEM_DIR}/sum_stats_final" \
  --pattern="${PATTERN}" \
  --ref="${GSEM_DIR}/geno/ref_gSEM_frq.txt" \
  --out_rdata="${OUT_RDATA}"

touch "${DONE_FILE}"
echo "[sumstats] Done."
