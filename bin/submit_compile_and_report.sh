#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"; shift 2;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash bin/submit_compile_and_report.sh --config config/my_run.sh" >&2
      exit 2;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "Usage: bash bin/submit_compile_and_report.sh --config config/my_run.sh" >&2
  exit 2
fi

CONFIG_FILE=$(realpath "${CONFIG_FILE}")
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: config not found: ${CONFIG_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

mkdir -p "${GSEM_DIR}/logs" || true

compile_jid=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --time="${COMPILE_TIME}" --cpus-per-task="${COMPILE_CPUS}" --mem="${COMPILE_MEM}" \
  --job-name="gsem_compile" \
  --output="${GSEM_DIR}/logs/%x_%j.out" \
  --error="${GSEM_DIR}/logs/%x_%j.err" \
  --chdir="${REPO_DIR}" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/compile_mlma.sbatch" | cut -d';' -f1)

report_jid=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --dependency="afterok:${compile_jid}" \
  --time="${REPORT_TIME}" --cpus-per-task="${REPORT_CPUS}" --mem="${REPORT_MEM}" \
  --job-name="gsem_report" \
  --output="${GSEM_DIR}/logs/%x_%j.out" \
  --error="${GSEM_DIR}/logs/%x_%j.err" \
  --chdir="${REPO_DIR}" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/prepare_report.sbatch" | cut -d';' -f1)

echo "Submitted compile/report chain:"
echo "  compile_mlma:  ${compile_jid}"
echo "  report:        ${report_jid}"
echo
echo "Monitor: squeue -u ${USER}"
echo "Logs: ${GSEM_DIR}/logs/<job>_<jobid>.out (.err)"
