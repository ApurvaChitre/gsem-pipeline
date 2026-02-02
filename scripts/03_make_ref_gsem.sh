#!/usr/bin/env bash
set -euo pipefail

# Generate reference allele frequencies for GenomicSEM sumstats()
#
# Usage:
#   bash scripts/03_make_ref_gsem.sh \
#     --bfile=/path/to/genotypes_prefix \
#     --out_prefix=/path/to/out/ref_gSEM \
#     --chr=1-20 \
#     --rscript=/path/to/scripts/03_make_ref_gSEM_frq_txt.R \
#     [--plink=/path/to/plink] \
#     [--force]

BFILE=""
OUT_PREFIX=""
CHR_RANGE="1-20"
RSCRIPT_PATH=""
PLINK_BIN="plink"
FORCE=0

for arg in "$@"; do
  case $arg in
    --bfile=*)
      BFILE="${arg#*=}"; shift;;
    --out_prefix=*)
      OUT_PREFIX="${arg#*=}"; shift;;
    --chr=*)
      CHR_RANGE="${arg#*=}"; shift;;
    --rscript=*)
      RSCRIPT_PATH="${arg#*=}"; shift;;
    --plink=*)
      PLINK_BIN="${arg#*=}"; shift;;
    --force)
      FORCE=1; shift;;
    *)
      # ignore unknown (allows passing through extra flags)
      ;;
  esac
done

if [[ -z "${BFILE}" || -z "${OUT_PREFIX}" || -z "${RSCRIPT_PATH}" ]]; then
  echo "ERROR: missing required args." >&2
  echo "  Required: --bfile= --out_prefix= --rscript=" >&2
  exit 2
fi

FRQ_FILE="${OUT_PREFIX}.frq"
FRQ_TXT="${OUT_PREFIX}_frq.txt"

mkdir -p "$(dirname "${OUT_PREFIX}")"

if [[ -f "${FRQ_TXT}" && ${FORCE} -eq 0 ]]; then
  echo "[ref] Found existing: ${FRQ_TXT} (use --force to regenerate)"
  exit 0
fi

if [[ ! -x "${PLINK_BIN}" ]]; then
  # allow PATH lookup
  if ! command -v "${PLINK_BIN}" >/dev/null 2>&1; then
    echo "ERROR: plink not found: ${PLINK_BIN}" >&2
    exit 2
  fi
fi

echo "[ref] Running plink --freq (chr=${CHR_RANGE})"
"${PLINK_BIN}" --bfile "${BFILE}" --chr "${CHR_RANGE}" --freq --out "${OUT_PREFIX}"

if [[ ! -f "${FRQ_FILE}" ]]; then
  echo "ERROR: Expected .frq file not found: ${FRQ_FILE}" >&2
  exit 1
fi

echo "[ref] Converting .frq -> $(basename "${FRQ_TXT}")"

# NOTE:
#   scripts/03_make_ref_gSEM_frq_txt.R expects:
#     --frq=/path/to/<prefix>.frq
#     --out=/path/to/<prefix>_frq.txt
#   (and optionally --bfile= or --bim= for SNP positions)
Rscript "${RSCRIPT_PATH}" \
  --frq="${FRQ_FILE}" \
  --bfile="${BFILE}" \
  --out="${FRQ_TXT}"

if [[ ! -s "${FRQ_TXT}" ]]; then
  echo "ERROR: ref freq output is missing/empty: ${FRQ_TXT}" >&2
  exit 1
fi

echo "[ref] Done: ${FRQ_TXT}"
