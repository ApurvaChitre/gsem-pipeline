#!/usr/bin/env bash
#SBATCH --job-name=gsem_prepare_report
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

mkdir -p "${MULTIVAR_REPORT_DIR}/results/gwas" \
         "${MULTIVAR_REPORT_DIR}/genome_info" \
         "${MULTIVAR_REPORT_DIR}/genotypes" \
         "${MULTIVAR_REPORT_DIR}/images" \
         "${MULTIVAR_REPORT_DIR}/log" \
         "${MULTIVAR_REPORT_DIR}/logerr" \
         "${MULTIVAR_REPORT_DIR}/temp" \
         "${MULTIVAR_REPORT_DIR}/data"

# 1) notebook template
if [[ -f "${REPO_DIR}/templates/multivariate_gwas_report/multivariate_gwas_report.ipynb" ]]; then
  cp -f "${REPO_DIR}/templates/multivariate_gwas_report/multivariate_gwas_report.ipynb" "${MULTIVAR_REPORT_DIR}/multivariate_gwas_report.ipynb"
fi

# 2) genome_info (optional)
if [[ -n "${GENOME_INFO_DIR:-}" && -d "${GENOME_INFO_DIR}" ]]; then
  rsync -av "${GENOME_INFO_DIR}/" "${MULTIVAR_REPORT_DIR}/genome_info/"
fi

# 3) genotypes directory: symlink the plink files from GENO_BFILE
#    (Some report pipelines expect genotypes/ to contain a plink prefix)
if [[ -n "${GENO_BFILE:-}" ]]; then
  GENO_DIR=$(dirname "${GENO_BFILE}")
  PREFIX=$(basename "${GENO_BFILE}")
  for ext in bed bim fam; do
    src="${GENO_DIR}/${PREFIX}.${ext}"
    dst="${MULTIVAR_REPORT_DIR}/genotypes/${PREFIX}.${ext}"
    if [[ -f "${src}" && ! -e "${dst}" ]]; then
      ln -s "${src}" "${dst}" || cp -f "${src}" "${dst}"
    fi
  done
fi

# 4) Optionally create the dummy processed_data_ready + data_dict files
if [[ "${RUN_CREATE_DUMMY_PHENO:-0}" == "1" ]]; then
  # Conda activation (recommended for R packages)
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

  Rscript "${REPO_DIR}/scripts/08_create_data_dict_processed_data_ready.R" \
    --measures="${REPORT_MEASURES}" \
    --processed_in="${PROCESSED_IN_FOR_REPORT}" \
    --data_dict_in="${DATA_DICT_IN_FOR_REPORT}" \
    --processed_out="${MULTIVAR_REPORT_DIR}/processed_data_ready.csv" \
    --data_dict_out="${MULTIVAR_REPORT_DIR}/data_dict_multivariate_gwas_report.csv" \
    --source_trait_col="${SOURCE_TRAIT_COL_FOR_REPORT}" \
    --source_dict_measure="${SOURCE_DICT_MEASURE_FOR_REPORT}"
fi

echo "[prepare_report] Done: ${MULTIVAR_REPORT_DIR}"
