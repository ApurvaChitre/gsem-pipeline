#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"; shift 2;;
    *)
      echo "Unknown argument: $1" >&2; exit 2;;
  esac
done

if [[ -z "${CONFIG_FILE}" ]]; then
  echo "Usage: bash bin/find_failed_usergwas_chunks.sh --config config/my_run.sh" >&2
  exit 2
fi

CONFIG_FILE=$(realpath "${CONFIG_FILE}")
source "${CONFIG_FILE}"

NUM_FILE="${GSEM_DIR}/split_sumstats/num_SNP_sets.txt"
if [[ ! -f "${NUM_FILE}" ]]; then
  echo "Missing: ${NUM_FILE} (run split_sumstats first)" >&2
  exit 1
fi

NUM=$(cat "${NUM_FILE}")
OUT_DIR="${GSEM_DIR}/results/multivariate_gwas/${MODEL_NAME}"

missing=()
for i in $(seq 1 "${NUM}"); do
  f="${OUT_DIR}/userGWAS_chunk_${i}.RData"
  if [[ ! -s "${f}" ]]; then
    missing+=("${i}")
  fi
done

if [[ ${#missing[@]} -eq 0 ]]; then
  echo "All chunks present (${NUM}/${NUM})."
  exit 0
fi

# Print as comma-separated list (Slurm --array accepts commas and ranges, but we keep it simple)
(IFS=,; echo "Missing chunks: ${missing[*]}")

