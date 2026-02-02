#!/usr/bin/env bash
#SBATCH --job-name=gsem_split_sumstats
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

mkdir -p "${GSEM_DIR}/split_sumstats" "${GSEM_DIR}/logs"

SUMSTATS_RDATA="${GSEM_DIR}/sum_stats_final/mySumstatsGSEM.RData"
DONE_FILE="${GSEM_DIR}/split_sumstats/.split.done"

if [[ "${SKIP_EXISTING}" == "1" && -f "${DONE_FILE}" && -f "${GSEM_DIR}/split_sumstats/num_SNP_sets.txt" ]]; then
  echo "[split] DONE_FILE exists; skipping."
  exit 0
fi

# Activate env (needs data.table)
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate "${CONDA_ENV}"

Rscript "${REPO_DIR}/scripts/05_split_sumstats.R" \
  --sumstats_rdata="${SUMSTATS_RDATA}" \
  --out_dir="${GSEM_DIR}/split_sumstats" \
  --subset_size="${SPLIT_SIZE}"

touch "${DONE_FILE}"
echo "[split] Done."
