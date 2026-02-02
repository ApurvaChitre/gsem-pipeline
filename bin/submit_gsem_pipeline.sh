#!/usr/bin/env bash

# MPH-style: fail fast on errors/pipe failures, but avoid `set -u` because
# (a) conda activation scripts sometimes reference unset vars
# (b) we prefer explicit required-variable checks for nicer messages
set -eo pipefail

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"; shift 2;;
    -h|--help)
      echo "Usage: bash bin/submit_gsem_pipeline.sh --config config/my_run.sh"; exit 0;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash bin/submit_gsem_pipeline.sh --config config/my_run.sh" >&2
      exit 2;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "Usage: bash bin/submit_gsem_pipeline.sh --config config/my_run.sh" >&2
  exit 2
fi

CONFIG_FILE=$(realpath "${CONFIG_FILE}")
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: config not found: ${CONFIG_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

# Required config
: "${GSEM_DIR:?Set GSEM_DIR in config}"
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT in config}"
: "${SLURM_PARTITION:?Set SLURM_PARTITION in config}"
: "${SLURM_QOS:?Set SLURM_QOS in config}"

mkdir -p "${GSEM_DIR}/logs" || true

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

# Convenience: only pass --time/--mem/--cpus-per-task if user set them in config.
# Otherwise Slurm will use the defaults baked into each *.sbatch file.
add_res_args() {
  local -n _arr=$1
  local _time="$2" _cpus="$3" _mem="$4"
  if [[ -n "${_time}" ]]; then _arr+=("--time=${_time}"); fi
  if [[ -n "${_cpus}" ]]; then _arr+=("--cpus-per-task=${_cpus}"); fi
  if [[ -n "${_mem}" ]]; then _arr+=("--mem=${_mem}"); fi
}

common=(
  --parsable
  --chdir="${REPO_DIR}"
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}"
  --output="${GSEM_DIR}/logs/%x_%j.out"
  --error="${GSEM_DIR}/logs/%x_%j.err"
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}"
)

# 0) Preflight: check MPH covstruc + quick 1-factor fit
# (Lets users catch non-positive definite matrices *before* spending time on sumstats/userGWAS.)
RUN_PREFLIGHT="${RUN_PREFLIGHT:-1}"
jid_preflight=""
dep_preflight=""
if [[ "${RUN_PREFLIGHT}" == "1" ]]; then
  args_pre=("${common[@]}" --job-name="gsem_preflight")
  add_res_args args_pre "${PREFLIGHT_TIME:-}" "${PREFLIGHT_CPUS:-}" "${PREFLIGHT_MEM:-}"
  jid_preflight=$(sbatch "${args_pre[@]}" "${REPO_DIR}/slurm/preflight_model.sbatch" | cut -d';' -f1)
  dep_preflight="afterok:${jid_preflight}"
fi

# 1) Make formatted sumstats from univariate GWAS
args_sum=("${common[@]}" --job-name="gsem_make_sumstats")
if [[ -n "${dep_preflight}" ]]; then
  args_sum+=(--dependency="${dep_preflight}")
fi
add_res_args args_sum "${MAKE_SUMSTATS_TIME:-}" "${MAKE_SUMSTATS_CPUS:-}" "${MAKE_SUMSTATS_MEM:-}"
jid_make_sumstats=$(sbatch "${args_sum[@]}" "${REPO_DIR}/slurm/make_sumstats.sbatch" | cut -d';' -f1)

# 2) Make reference freq file for GenomicSEM
args_ref=("${common[@]}" --job-name="gsem_make_ref" --dependency="afterok:${jid_make_sumstats}")
add_res_args args_ref "${MAKE_REF_TIME:-}" "${MAKE_REF_CPUS:-}" "${MAKE_REF_MEM:-}"
jid_make_ref=$(sbatch "${args_ref[@]}" "${REPO_DIR}/slurm/make_ref.sbatch" | cut -d';' -f1)

# 3) Run sumstats() and write mySumstatsGSEM.RData
args_ss=("${common[@]}" --job-name="gsem_sumstats" --dependency="afterok:${jid_make_ref}")
add_res_args args_ss "${SUMSTATS_TIME:-}" "${SUMSTATS_CPUS:-}" "${SUMSTATS_MEM:-}"
jid_sumstats=$(sbatch "${args_ss[@]}" "${REPO_DIR}/slurm/sumstats.sbatch" | cut -d';' -f1)

# 4) Split sumstats RData into chunks for userGWAS
args_split=("${common[@]}" --job-name="gsem_split_sumstats" --dependency="afterok:${jid_sumstats}")
add_res_args args_split "${SPLIT_TIME:-}" "${SPLIT_CPUS:-}" "${SPLIT_MEM:-}"
jid_split=$(sbatch "${args_split[@]}" "${REPO_DIR}/slurm/split_sumstats.sbatch" | cut -d';' -f1)

# 5) Launch Phase A (usermodel fit) and optionally Phase B (userGWAS + report)
# NOTE: This job *runs the model fit* inside it, so give it enough memory/time.
args_launch=("${common[@]}" --job-name="gsem_launch" --dependency="afterok:${jid_split}")
# Prefer USERMODEL_* resources for the launch job, since it runs the fit.
add_res_args args_launch "${USERMODEL_TIME:-}" "${USERMODEL_CPUS:-}" "${USERMODEL_MEM:-}"
jid_launch=$(sbatch "${args_launch[@]}" "${REPO_DIR}/slurm/launch_usergwas_chain.sbatch" | cut -d';' -f1)

cat <<MSG
Submitted Genomic SEM pipeline.

Job IDs:
  preflight    : ${jid_preflight}
  make_sumstats : ${jid_make_sumstats}
  make_ref      : ${jid_make_ref}
  sumstats      : ${jid_sumstats}
  split_sumstats: ${jid_split}
  launch        : ${jid_launch}

Logs:
  ${GSEM_DIR}/logs/

Phase A outputs (after launch finishes):
  ${GSEM_DIR}/results/model_fit/

If RUN_USERGWAS=0 in your config, the pipeline will STOP after model fit.
To proceed to Phase B (userGWAS + report), set RUN_USERGWAS=1 and re-run:
  bash bin/submit_gsem_pipeline.sh --config ${CONFIG_FILE}

(Important: when copying commands, use normal hyphens "--config" (two dashes), not an en-dash.)
MSG
