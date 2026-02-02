#!/usr/bin/env bash
set -eo pipefail

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"; shift 2;;
    -h|--help)
      echo "Usage: bash bin/submit_debug_env.sh --config config/my_run.sh"; exit 0;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash bin/submit_debug_env.sh --config config/my_run.sh" >&2
      exit 2;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "Usage: bash bin/submit_debug_env.sh --config config/my_run.sh" >&2
  exit 2
fi

CONFIG_FILE=$(realpath "${CONFIG_FILE}")
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: config not found: ${CONFIG_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"
: "${GSEM_DIR:?Set GSEM_DIR in config}"

mkdir -p "${GSEM_DIR}/logs" || true

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

jid=$(sbatch --parsable \
  --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" --qos="${SLURM_QOS}" \
  --time="0:15:00" --cpus-per-task=1 --mem=4G \
  --chdir="${REPO_DIR}" \
  --job-name="gsem_debug_env" \
  --output="${GSEM_DIR}/logs/%x_%j.out" \
  --error="${GSEM_DIR}/logs/%x_%j.err" \
  --export=ALL,CONFIG_FILE="${CONFIG_FILE}",REPO_DIR="${REPO_DIR}" \
  "${REPO_DIR}/slurm/debug_env.sbatch" | cut -d';' -f1)

echo "Submitted debug job: ${jid}"
echo "Logs: ${GSEM_DIR}/logs/gsem_debug_env_${jid}.out (.err)"
