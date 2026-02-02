#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE=""
ARRAY_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"; shift 2;;
    --array)
      ARRAY_OVERRIDE="$2"; shift 2;;
    --array=*)
      ARRAY_OVERRIDE="${1#*=}"; shift 1;;
    -h|--help)
      cat <<EOF
Usage:
  bash bin/submit_usergwas_array.sh --config config/my_run.sh [--array "1-100%20"|"5,9,12"]

Notes:
  - If --array is omitted, we read num_SNP_sets.txt and submit 1-N (optionally throttled by USERGWAS_ARRAY_MAX).
  - If --array is provided, it is passed directly to sbatch --array (useful for rerunning failed chunks).
EOF
      exit 0;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash bin/submit_usergwas_array.sh --config config/my_run.sh [--array \"...\"]" >&2
      exit 2;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "Usage: bash bin/submit_usergwas_array.sh --config config/my_run.sh [--array \"...\"]" >&2
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

# Determine array spec
if [[ -n "${ARRAY_OVERRIDE}" ]]; then
  ARRAY_SPEC="${ARRAY_OVERRIDE}"
else
  NUM_FILE="${GSEM_DIR}/split_sumstats/num_SNP_sets.txt"
  if [[ ! -f "${NUM_FILE}" ]]; then
    echo "ERROR: missing ${NUM_FILE}. Run split_sumstats first." >&2
    exit 2
  fi

  NUM=$(cat "${NUM_FILE}")
  if ! [[ "${NUM}" =~ ^[0-9]+$ ]] || [[ "${NUM}" -lt 1 ]]; then
    echo "ERROR: invalid chunk count in ${NUM_FILE}: ${NUM}" >&2
    exit 2
  fi

  # Optional throttle: max number of concurrent array tasks
  ARRAY_MAX="${USERGWAS_ARRAY_MAX:-}"
  ARRAY_SPEC="1-${NUM}"
  if [[ -n "${ARRAY_MAX}" ]]; then
    if [[ "${ARRAY_MAX}" =~ ^[0-9]+$ ]] && [[ "${ARRAY_MAX}" -gt 0 ]]; then
      ARRAY_SPEC="1-${NUM}%${ARRAY_MAX}"
    fi
  fi
fi

jid=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --time="${USERGWAS_TIME}" --cpus-per-task="${USERGWAS_CPUS}" --mem="${USERGWAS_MEM}" \
  --array="${ARRAY_SPEC}" \
  --job-name="gsem_usergwas" \
  --chdir="${REPO_DIR}" \
  --output="${GSEM_DIR}/logs/%x_%A_%a.out" \
  --error="${GSEM_DIR}/logs/%x_%A_%a.err" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/usergwas_array.sbatch" | cut -d';' -f1)

echo "Submitted userGWAS array: ${jid}"
echo "Array spec: ${ARRAY_SPEC}"
echo "Logs: ${GSEM_DIR}/logs/gsem_usergwas_${jid}_<TASKID>.out (.err)"
