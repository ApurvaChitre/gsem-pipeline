#!/usr/bin/env bash
#SBATCH --job-name=gsem_launch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00

# Launch job: runs the usermodel fit (Phase A) and then submits Phase B jobs (optional).
# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

mkdir -p "${GSEM_DIR}/logs" "${GSEM_DIR}/results/model_fit"

# Activate conda env (R + packages)
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate "${CONDA_ENV}"

echo "[env] Activated conda env: ${CONDA_ENV}"
echo "[env] Rscript: $(command -v Rscript)"
echo "[env] plink:   $(command -v plink || true)"

# ------------------------------
# Phase A: usermodel fit (recommended)
# ------------------------------
if [[ "${RUN_USERMODEL_FIT}" == "1" ]]; then
  echo "[launch] Phase A: fitting usermodel: ${MODEL_NAME}"

  MODEL_ARGS=(--model_name "${MODEL_NAME}" --out_dir "${GSEM_DIR}/results/model_fit")
  if [[ -n "${MODEL_FILE:-}" ]]; then
    MODEL_ARGS+=(--model_file "${MODEL_FILE}")
  fi

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
    "${DICT_ARG[@]}" \
    "${MODEL_ARGS[@]}"

  echo "[launch] Model fit outputs: ${GSEM_DIR}/results/model_fit/"
else
  echo "[launch] RUN_USERMODEL_FIT=0; skipping model fit."
fi

# Stop here unless user explicitly wants userGWAS.
if [[ "${RUN_USERGWAS}" != "1" ]]; then
  cat <<MSG
[launch] RUN_USERGWAS != 1
[launch] STOPPING after Phase A.

Next steps:
  1) Inspect model fit outputs:
     ${GSEM_DIR}/results/model_fit/

  2) If the model looks good and you want to run userGWAS, set RUN_USERGWAS=1 in your config
     and re-submit:
       bash bin/submit_gsem_pipeline.sh --config ${CONFIG_FILE}
MSG
  exit 0
fi

# ------------------------------
# Phase B: userGWAS array + compile + report
# ------------------------------
mkdir -p "${GSEM_DIR}/results/multivariate_gwas" "${GSEM_DIR}/multivariate_gwas_report"

# Submit userGWAS array
EXPECTED_FILE="${GSEM_DIR}/split_sumstats/num_SNP_sets.txt"
if [[ ! -f "${EXPECTED_FILE}" ]]; then
  echo "ERROR: expected chunk count file not found: ${EXPECTED_FILE}" >&2
  echo "Did split_sumstats finish successfully?" >&2
  exit 2
fi

SPLIT_NCHUNKS=$(cat "${EXPECTED_FILE}")
if ! [[ "${SPLIT_NCHUNKS}" =~ ^[0-9]+$ ]] || [[ "${SPLIT_NCHUNKS}" -lt 1 ]]; then
  echo "ERROR: invalid chunk count in ${EXPECTED_FILE}: ${SPLIT_NCHUNKS}" >&2
  exit 2
fi

# Optional throttle: max number of concurrent array tasks
ARRAY_MAX="${USERGWAS_ARRAY_MAX:-}"
ARRAY_SPEC="1-${SPLIT_NCHUNKS}"
if [[ -n "${ARRAY_MAX}" ]]; then
  if [[ "${ARRAY_MAX}" =~ ^[0-9]+$ ]] && [[ "${ARRAY_MAX}" -gt 0 ]]; then
    ARRAY_SPEC="1-${SPLIT_NCHUNKS}%${ARRAY_MAX}"
  fi
fi

jid_array=$(sbatch --parsable   --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}"   --time="${USERGWAS_TIME}" --cpus-per-task="${USERGWAS_CPUS}" --mem="${USERGWAS_MEM}"   --array="${ARRAY_SPEC}"   --job-name="gsem_usergwas"   --output="${GSEM_DIR}/logs/%x_%A_%a.out"   --error="${GSEM_DIR}/logs/%x_%A_%a.err"   --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}"   "${REPO_DIR}/slurm/usergwas_array.sbatch" | cut -d';' -f1)

echo "[launch] Submitted userGWAS array: ${jid_array}"

# Compile userGWAS outputs into MLMA-like format
jid_compile=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --time="${COMPILE_TIME}" --cpus-per-task="${COMPILE_CPUS}" --mem="${COMPILE_MEM}" \
  --dependency="afterok:${jid_array}" \
  --job-name="gsem_compile" \
  --output="${GSEM_DIR}/logs/%x_%j.out" \
  --error="${GSEM_DIR}/logs/%x_%j.err" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/compile_mlma.sbatch" | cut -d';' -f1)

echo "[launch] Submitted compile step: ${jid_compile}"

# Prepare multivariate_gwas_report folder scaffold + helper files
jid_report=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --time="${REPORT_TIME}" --cpus-per-task="${REPORT_CPUS}" --mem="${REPORT_MEM}" \
  --dependency="afterok:${jid_compile}" \
  --job-name="gsem_report" \
  --output="${GSEM_DIR}/logs/%x_%j.out" \
  --error="${GSEM_DIR}/logs/%x_%j.err" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/prepare_report.sbatch" | cut -d';' -f1)

echo "[launch] Submitted report-prep step: ${jid_report}"

echo "[launch] Done submitting Phase B jobs."
