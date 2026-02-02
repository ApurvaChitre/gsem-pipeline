#!/usr/bin/env bash
#SBATCH --job-name=gsem_make_ref
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

# ----- Conda activation (recommended: provides plink + R packages) -----
if [[ -n "${CONDA_ENV:-}" ]]; then
  if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
  else
    [[ -f "${HOME}/.bashrc" ]] && source "${HOME}/.bashrc"
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
  fi
  conda activate "${CONDA_ENV}"
  echo "[env] Activated conda env: ${CONDA_ENV}"
  echo "[env] Rscript: $(command -v Rscript || true)"
  echo "[env] plink:   $(command -v plink || true)"
fi


mkdir -p "${GSEM_DIR}/geno" "${GSEM_DIR}/logs"

REF_TXT="${GSEM_DIR}/geno/ref_gSEM_frq.txt"
DONE_FILE="${GSEM_DIR}/geno/.make_ref.done"

if [[ "${SKIP_EXISTING}" == "1" && -f "${DONE_FILE}" && -s "${REF_TXT}" ]]; then
  echo "[make_ref] DONE_FILE exists and ref file non-empty; skipping."
  exit 0
fi

echo "[make_ref] Writing ref allele freqs to: ${REF_TXT}"

bash "${REPO_DIR}/scripts/03_make_ref_gsem.sh" \
  --bfile="${GENO_BFILE}" \
  --out_prefix="${GSEM_DIR}/geno/ref_gSEM" \
  --chr="${CHR_RANGE}" \
  --plink="${PLINK_BIN}" \
  --rscript="${REPO_DIR}/scripts/03_make_ref_gSEM_frq_txt.R"

touch "${DONE_FILE}"
echo "[make_ref] Done."
