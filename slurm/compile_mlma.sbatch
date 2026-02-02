#!/usr/bin/env bash
#SBATCH --job-name=gsem_compile_mlma
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

IN_DIR="${GSEM_DIR}/results/multivariate_gwas/${MODEL_NAME}"
OUT_DIR="${MULTIVAR_REPORT_DIR}/results/gwas"
mkdir -p "${OUT_DIR}"

# Safety: refuse to compile if any expected chunk output is missing
EXPECTED_FILE="${GSEM_DIR}/split_sumstats/num_SNP_sets.txt"
if [[ -f "${EXPECTED_FILE}" ]]; then
  expected=$(cat "${EXPECTED_FILE}")
  missing=0
  for i in $(seq 1 "${expected}"); do
    f="${IN_DIR}/userGWAS_chunk_${i}.RData"
    if [[ ! -s "${f}" ]]; then
      echo "MISSING chunk output: ${f}" >&2
      missing=1
    fi
  done
  if [[ "${missing}" -eq 1 ]]; then
    echo "Refusing to compile: missing chunk outputs. Resubmit failed chunks first." >&2
    exit 2
  fi
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

Rscript "${REPO_DIR}/scripts/07_compile_rdata_to_mlma.R" \
  --input_dir="${IN_DIR}" \
  --out_dir="${OUT_DIR}" \
  --prefix="regressedlr_${MODEL_NAME}"
