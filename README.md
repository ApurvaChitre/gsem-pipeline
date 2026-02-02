# Genomic SEM usermodel (1-factor) pipeline

This repo is a **TSCC/Slurm-friendly** workflow to run a **GenomicSEM usermodel GWAS with a single common factor (1F)**, then convert outputs into the folder layout expected by the downstream `multivariate_gwas_report.ipynb` workflow.

It is intentionally narrow:
- ✅ supports **one-factor usermodel** (common factor)
- ✅ handles **any number of traits** (no hard-coded “N=6 traits” assumptions)
- ✅ creates **`results/gwas/*.mlma`** outputs for downstream reporting
- ✅ includes helpers to **resubmit only failed GWAS chunks**
- ❌ does **not** ship code for “1F + residual covariances”, 2-factor variants, bifactor, clustering/heatmaps, etc.

---

## What you need before running

1. **MPH pipeline outputs** (from the separate MPH pipeline repo):
   - `${MPH_DIR}/gsem/MPH_genomicSEM.RData`
   - `${MPH_DIR}/pheno/pheno_cohort_project_dict.csv`

2. **Univariate GWAS results** for each trait (same naming as MPH traits, typically `regressedlr_*`).
   - You point the pipeline to the directory that contains `*.mlma` files (recursively).

3. A working **GenomicSEM conda environment** (you can reuse the same env as MPH).

---

## Quick start (TSCC)

### 1) Clone repo

```
git clone https://github.com/<YOU>/<REPO>.git
cd <REPO>
```

### 2) Conda environment

If you already made the `genomic_sem` env while running the MPH pipeline, just:

```
conda activate genomic_sem
```

Otherwise:
- Use the **lock file from the MPH pipeline repo** (recommended for reproducibility), OR
- Use this repo's lightweight example file.

If you have **mamba**:

```
mamba env create -f env/environment.example.yml -n genomic_sem
```

If you only have **conda**:

```
conda env create -f env/environment.example.yml -n genomic_sem
```

If `library(GenomicSEM)` fails after that, run:

```
Rscript env/install_genomicsem_from_github.R
```

### 3) Create your config

Copy and edit:

```
cp config/example_run.sh config/my_run.sh
nano config/my_run.sh
```

### 4) Debug env + paths (strongly recommended)

This runs a short Slurm job that verifies:
- required MPH outputs exist
- your univariate GWAS directory has `*.mlma`
- `plink`, `Rscript`, and `GenomicSEM` are available after conda activation

```
bash bin/submit_debug_env.sh --config config/my_run.sh

Note: if you copy/paste from Slack/Google Docs and see `Unknown argument: –config`,
that’s an **en-dash**. Re-type it as **two normal hyphens**: `--config`.
```

### 5) Submit the full pipeline

This repo is designed to prevent accidental waste: **the userGWAS array is the expensive step**.

Recommended workflow is:

**Step 5A: run through model-fit only (safe / cheap)**

1) In your config, keep:

RUN_USERGWAS="0"

2) Submit:

bash bin/submit_gsem_pipeline.sh --config config/my_run.sh

This will:
- format sumstats
- split sumstats into SNP chunks
- fit the 1-factor usermodel and write outputs under `${GSEM_DIR}/results/model_fit/`
- **STOP** (no userGWAS array launched)

**Step 5B: if model fit looks OK, launch the full userGWAS array**

1) In your config, set:

RUN_USERGWAS="1"

2) Re-submit (completed steps will be skipped):

bash bin/submit_gsem_pipeline.sh --config config/my_run.sh

Monitor:

```
squeue -u $USER
```

---

## Outputs

Key outputs (paths are all under `${GSEM_DIR}`):

- `sum_stats_final/` : formatted per-trait sumstats + `mySumstatsGSEM.RData`
- `geno/` : reference allele frequency file (`ref_gSEM_frq.txt`)
- `split_sumstats/subsets/` : chunked SNP sets (`sumstats_subset_*.RData`)
- `results/multivariate_gwas/${MODEL_NAME}/` : per-chunk `chunk_*.RData` userGWAS outputs
- `multivariate_gwas_report/results/gwas/` : compiled `*.mlma` outputs (ready for downstream reporting)

---

## Resuming after failures (common in GenomicSEM)

If the userGWAS array stops due to a chunk failing (e.g., Heywood / convergence issues), you typically want to **rerun only the missing chunks**:

1) Find missing chunks:

```
bash bin/find_missing_usergwas_chunks.sh --config config/my_run.sh
```

2) Resubmit only those indices:

```
bash bin/submit_usergwas_array.sh --config config/my_run.sh --array "12,45,78"
```

3) When all chunks exist, run compile + report steps:

```
bash bin/submit_compile_and_report.sh --config config/my_run.sh
```

---

## Notes on model specification

This pipeline runs a **single 1-factor usermodel**.

By default it auto-builds:

```
F1 =~ trait1 + trait2 + ... + traitN
```

You *can* override the model string by pointing `MODEL_FILE` to a text file containing a lavaan-style model. This is mainly intended for minor tweaks; it’s still on you to validate the model fit.

---

## Repo layout

- `config/` : example config + (optional) model text file
- `bin/` : submission + helper scripts
- `slurm/` : sbatch scripts
- `scripts/` : R + bash pipeline logic
- `templates/multivariate_gwas_report/` : notebook template + expected folder skeleton

