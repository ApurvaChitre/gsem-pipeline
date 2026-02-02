#!/usr/bin/env bash
#SBATCH --job-name=gsem_usergwas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is not set" >&2
  exit 2
fi

# Inputs
MPH_RDATA="${MPH_DIR}/gsem/MPH_genomicSEM.RData"

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
SUBSET_RDATA="${GSEM_DIR}/split_sumstats/sumstats_subset_${SLURM_ARRAY_TASK_ID}.RData"

# Outputs
OUT_DIR="${GSEM_DIR}/results/multivariate_gwas/${MODEL_NAME}"
mkdir -p "${OUT_DIR}"
OUT_RDATA="${OUT_DIR}/userGWAS_chunk_${SLURM_ARRAY_TASK_ID}.RData"

if [[ "${SKIP_EXISTING}" == "1" && -s "${OUT_RDATA}" ]]; then
  echo "[skip] ${OUT_RDATA} exists"
  exit 0
fi

if [[ ! -s "${SUBSET_RDATA}" ]]; then
  echo "Missing subset RData: ${SUBSET_RDATA}" >&2
  exit 1
fi

# Conda activation (optional)
if [[ -n "${CONDA_ENV}" ]]; then
  if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
    source "${CONDA_SH}"
  else
    BASE=$(conda info --base 2>/dev/null || true)
    if [[ -n "${BASE}" && -f "${BASE}/etc/profile.d/conda.sh" ]]; then
      source "${BASE}/etc/profile.d/conda.sh"
    fi
  fi
  conda activate "${CONDA_ENV}"
fi

Rscript "${REPO_DIR}/scripts/06_usergwas_chunk.R" \
  --mph_rdata "${MPH_RDATA}" \
  --dict_file "${DICT_F}" \
  --subset_rdata "${SUBSET_RDATA}" \
  --out_rdata "${OUT_RDATA}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --model_file "${MODEL_FILE}"
