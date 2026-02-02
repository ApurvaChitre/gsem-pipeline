#!/usr/bin/env bash
#SBATCH --job-name=gsem_debug_env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:10:00

# "Soft strict" mode: fail on errors/pipe failures, but do NOT die on unset vars.
# Conda activate/deactivate scripts may reference unset vars (e.g. CONDA_BACKUP_CXX)
# and will crash under `set -u`.
set -eo pipefail

: "${CONFIG_FILE?CONFIG_FILE not set (pass via --export)}"
source "${CONFIG_FILE}"

# Where is the repo? (this file lives in repo/slurm/)
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_DIR}"

# Logging helper
fail=0
say() { echo "$1"; }

say "=== GSEM DEBUG ENV ==="
say "host=$(hostname)"
say "date=$(date)"
say "pwd=$(pwd)"
say "SLURM_JOB_ID=${SLURM_JOB_ID:-<none>}"
say "CONFIG_FILE=${CONFIG_FILE}"

say ""
say "== Config summary =="
say "MPH_DIR=${MPH_DIR}"
say "GSEM_DIR=${GSEM_DIR}"
say "UNIV_GWAS_DIR=${UNIV_GWAS_DIR}"
say "GENO_BFILE=${GENO_BFILE}"
say "CHR_RANGE=${CHR_RANGE}"
say "CONDA_ENV=${CONDA_ENV:-<unset>}"
say "CONDA_SH=${CONDA_SH:-<auto>}"

say ""
say "== Basic path checks =="
# Outputs
if mkdir -p "${GSEM_DIR}/logs" 2>/dev/null; then
  say "OK: can write GSEM_DIR -> ${GSEM_DIR}"
else
  say "FAIL: cannot create GSEM_DIR -> ${GSEM_DIR}"; fail=1
fi

# MPH inputs
mph_rdata="${MPH_DIR}/gsem/MPH_genomicSEM.RData"

# Dict file path (trait list + N). Prefer explicit DICT_FILE in config; fallback to common MPH layouts.
dict_file="${DICT_FILE:-}"
if [[ -z "${dict_file}" ]]; then
  dict_file="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"
fi
if [[ ! -f "${dict_file}" ]]; then
  for cand in \
    "${MPH_DIR}/pheno/pheno_cohort_project_dict.csv" \
    "${MPH_DIR}/pheno_cohort_project_dict.csv" \
    "${MPH_DIR}/individual_projects_seq/pheno_cohort_project_dict.csv" \
    ; do
    if [[ -f "${cand}" ]]; then
      dict_file="${cand}"
      break
    fi
  done
fi
if [[ -f "${mph_rdata}" ]]; then say "OK: MPH_genomicSEM -> ${mph_rdata}"; else say "FAIL: missing MPH_genomicSEM -> ${mph_rdata}"; fail=1; fi
if [[ -f "${dict_file}" ]]; then say "OK: pheno dict -> ${dict_file}"; else say "FAIL: missing pheno dict -> ${dict_file}"; fail=1; fi

# Univariate GWAS dir
if [[ -d "${UNIV_GWAS_DIR}" ]]; then
  say "OK: UNIV_GWAS_DIR exists"
else
  say "FAIL: UNIV_GWAS_DIR missing -> ${UNIV_GWAS_DIR}"; fail=1
fi

# Genotypes
if [[ -f "${GENO_BFILE}.bed" && -f "${GENO_BFILE}.bim" && -f "${GENO_BFILE}.fam" ]]; then
  say "OK: GENO_BFILE bed/bim/fam exist"
else
  say "FAIL: GENO_BFILE missing bed/bim/fam -> ${GENO_BFILE}.*"; fail=1
fi

say ""
say "== Tooling checks (before conda) =="
for cmd in awk cut wc perl; do
  if command -v "$cmd" >/dev/null 2>&1; then
    say "OK: $cmd -> $(command -v $cmd)"
  else
    say "FAIL: missing command '$cmd'"; fail=1
  fi
done

# plink may be in conda env
if command -v plink >/dev/null 2>&1; then
  say "OK: plink -> $(command -v plink)"
else
  say "WARN: plink not on PATH (expected if you rely on conda env)"
fi

if command -v Rscript >/dev/null 2>&1; then
  say "OK: Rscript -> $(command -v Rscript)"
else
  say "WARN: Rscript not on PATH (expected if you rely on conda env)"
fi

say ""
say "== Conda activation (optional) =="
if [[ -n "${CONDA_ENV:-}" ]]; then
  # Find conda.sh
  conda_sh="${CONDA_SH:-}"
  if [[ -z "${conda_sh}" ]]; then
    # best-effort: conda might be on PATH
    if command -v conda >/dev/null 2>&1; then
      conda_base=$(conda info --base)
      conda_sh="${conda_base}/etc/profile.d/conda.sh"
    fi
  fi

  if [[ -n "${conda_sh}" && -f "${conda_sh}" ]]; then
    # shellcheck source=/dev/null
    source "${conda_sh}"
    conda activate "${CONDA_ENV}"
    say "OK: activated conda env: ${CONDA_ENV}"
    say "Rscript -> $(command -v Rscript 2>/dev/null || echo '<missing>')"
    say "plink   -> $(command -v plink   2>/dev/null || echo '<missing>')"

    # R package sanity
    Rscript - <<'RS'
pkgs <- c('GenomicSEM','data.table','readr','dplyr','stringr','optparse')
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(missing)) {
  cat('MISSING R packages:', paste(missing, collapse=', '), '\n')
  quit(status=2)
}
cat('OK: required R packages are installed\n')
RS
    say "OK: required R packages are installed"
  else
    say "FAIL: could not locate conda.sh. Set CONDA_SH in your config."; fail=1
  fi
else
  say "NOTE: CONDA_ENV is empty; skipping conda activation checks"
fi

say ""
if [[ "${fail}" -eq 0 ]]; then
  say "DEBUG RESULT: PASS"
  say "=== DONE ==="
  exit 0
else
  say "DEBUG RESULT: FAIL"
  say "=== DONE ==="
  exit 1
fi
