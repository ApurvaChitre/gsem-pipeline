# Genomic SEM usermodel (1-factor) pipeline

This repo is a **TSCC/Slurm-friendly** workflow to run a **GenomicSEM usermodel GWAS** (default: a single common factor / 1-factor model), then convert outputs into the folder layout expected by the downstream `multivariate_gwas_report.ipynb` workflow.

It is intentionally narrow:

- ships a default **1-factor usermodel** (common factor), but supports a **custom user-supplied model** via `MODEL_FILE` (including multi-factor)
- handles **any number of traits** (no hard-coded “N=6 traits” assumptions)
- creates **`results/gwas/*.mlma`** outputs for downstream reporting
- includes helpers to **resubmit only failed GWAS chunks**
- does not include a library of pre-built residual-covariance / bifactor / multi-factor templates; custom models are supported but you must supply and validate them

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
- **stop** (no userGWAS array launched)

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

## How your model names become MLMA filenames and report measures

The pipeline compiles userGWAS chunk `.RData` outputs into **one `.mlma` per chromosome per target** (target = whatever is regressed on `SNP`).

### 1) Where the `<target>` / `<factor>` name comes from

For **usermodel/userGWAS** outputs, the compiler looks for rows where:

- `op == "~"`
- `rhs == "SNP"`

Every line in your model like:

```
<TARGET> ~ SNP
```

produces an MLMA output set for `<TARGET>`.

So:

- `F1 ~ SNP` → outputs labeled `f1`
- `F2 ~ SNP` → outputs labeled `f2`
- `Internalizing ~ SNP` → outputs labeled `internalizing`
- `addiction_1 ~ SNP` → outputs labeled `addiction_1`

If your custom model has **two latent factors** but you only include `F1 ~ SNP` (and not `F2 ~ SNP`), you will only get outputs for `F1`.

### 2) How the compiler converts names into filenames

The compiled MLMA files are written to:

`{MULTIVAR_REPORT_DIR}/results/gwas/`

Basename pattern:

`regressedlr_<MODEL_NAME>_<modelgroup>_<target>_chrgwas<CHR>.mlma`

Notes:

- `<modelgroup>` is the chunk-file group name. For this pipeline’s userGWAS chunk outputs it is `userGWAS`, which becomes `usergwas` in the filename.
- The MLMA basename is forced to **lowercase**.
- `<target>` is “slugged” for filenames: any character outside `[A-Za-z0-9._-]` is replaced with `_`, then lowercased.

Practical advice: **use simple factor names** like `F1`, `F2`, `internalizing`, `addiction_1`. Avoid spaces and punctuation in model labels.

### 3) Examples

#### Example A: default 1-factor model

Model contains:

```
F1 =~ trait1 + trait2 + trait3
F1 ~ SNP
```

Outputs include:

- `regressedlr_cf_resid_both_usergwas_f1_chrgwas1.mlma`
- ...
- `regressedlr_cf_resid_both_usergwas_f1_chrgwas20.mlma`

Trait string used by the report notebook (basename without `_chrgwasX.mlma`):

- `regressedlr_cf_resid_both_usergwas_f1`

Matching `data_dict` measure (remove `regressedlr_` prefix):

- `cf_resid_both_usergwas_f1`

#### Example B: 2-factor model (F1 and F2)

Model contains:

```
F1 =~ trait1 + trait2 + trait3
F2 =~ trait4 + trait5 + trait6
F1 ~ SNP
F2 ~ SNP
```

You will get two sets of outputs:

- `regressedlr_cf_resid_both_usergwas_f1_chrgwas*.mlma`
- `regressedlr_cf_resid_both_usergwas_f2_chrgwas*.mlma`

So you have **two measures**:

- `cf_resid_both_usergwas_f1`
- `cf_resid_both_usergwas_f2`

#### Example C: custom target names

Model contains:

```
Internalizing =~ trait1 + trait2 + trait3
Externalizing =~ trait4 + trait5 + trait6
Internalizing ~ SNP
Externalizing ~ SNP
```

Outputs will include:

- `regressedlr_cf_resid_both_usergwas_internalizing_chrgwas*.mlma`
- `regressedlr_cf_resid_both_usergwas_externalizing_chrgwas*.mlma`

Measures:

- `cf_resid_both_usergwas_internalizing`
- `cf_resid_both_usergwas_externalizing`

### 4) What to put in processed_data_ready.csv vs data_dict vs REPORT_MEASURES

These three MUST match the compiled naming:

- **MLMA trait string** (what the notebook keys on):
  - `regressedlr_<MODEL_NAME>_usergwas_<target>`
- **processed_data_ready.csv column name**:
  - same as the trait string above
- **data_dict measure**:
  - same string but with the leading `regressedlr_` removed
- **REPORT_MEASURES in config**:
  - comma-separated list of those data_dict measures (no `regressedlr_`)

Example (2-factor F1/F2):

```
REPORT_MEASURES="cf_resid_both_usergwas_f1,cf_resid_both_usergwas_f2"
```

### 5) If you’re unsure: auto-extract the measures from the compiled MLMA files

After compile finishes:

```
MLMA_DIR="${MULTIVAR_REPORT_DIR}/results/gwas"

ls "${MLMA_DIR}"/*.mlma \
  | sed -E 's#.*/##' \
  | sed -E 's/^regressedlr_//' \
  | sed -E 's/_chrgwas[0-9]+\.mlma$//' \
  | sort -u
```

Those printed strings are exactly what should go into `REPORT_MEASURES` (comma-join them).

---

## Creating a dummy processed_data_ready.csv and data_dict for the report (optional)

Some report workflows expect two small metadata files in the project folder:

- `processed_data_ready.csv` (must contain the trait column, i.e. `regressedlr_<...>`)
- a `data_dict_*.csv` with a `measure` column (typically *without* the `regressedlr_` prefix)

You have two options.

### Option A: Have the pipeline generate them (copy-from-mega approach)

If you have a “mega” `processed_data_ready.csv` and `data_dict_mega.csv`, you can copy a single existing trait column into one or more new dummy measures.

In your config (example):

```
RUN_CREATE_DUMMY_PHENO="1"

# Comma-separated list. You may include or omit the leading regressedlr_.
# (The script normalizes to data_dict-style measures without the prefix.)
REPORT_MEASURES="cf_resid_both_usergwas_f1,cf_resid_both_usergwas_f2"

# Inputs to copy from
PROCESSED_IN="/path/to/mega/processed_data_ready.csv"
DATA_DICT_IN="/path/to/mega/data_dict_mega.csv"
SOURCE_TRAIT_COL="regressedlr_locomotor_mega"
SOURCE_DICT_MEASURE="locomotor_mega"
```

When you run the compile/report step, it will write:

- `${MULTIVAR_REPORT_DIR}/processed_data_ready.csv`
- `${MULTIVAR_REPORT_DIR}/data_dict_multivariate_gwas_report.csv`

### Option B: Let the notebook create minimal files

The notebook template includes a small block that creates minimal versions from the genotype `.fam` file if they are missing. Trait values are dummy zeros (sufficient for report scaffolding).

---

## Example: plotting a 1-factor path diagram

See `examples/path_diagram_1factor_example.R` for a minimal example that:

- loads `MPH_genomicSEM.RData`
- builds a 1-factor model string
- fits `GenomicSEM::usermodel()`
- renders a path diagram with `semPlot`

---

## Resuming after failures (common in GenomicSEM)

If the userGWAS array stops due to a chunk failing (e.g., Heywood / convergence issues), you typically want to rerun only the missing chunks:

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

This pipeline ships a **default 1-factor usermodel**.

By default it auto-builds:

```
F1 =~ trait1 + trait2 + ... + traitN
F1 ~ SNP
```

You can override the model by pointing `MODEL_FILE` to a text file containing a lavaan-style model.

- Your custom model may be **multi-factor** (e.g., `F1`, `F2`, …).
- You will only get MLMA outputs for targets that have an explicit `~ SNP` line (e.g., `F2 ~ SNP`).
- The `<target>` name you use in the model (left-hand side of `~ SNP`) propagates into the compiled MLMA filenames (after slugging + lowercasing). See the section above on naming/`REPORT_MEASURES`.

Custom models are powerful, but it’s still on you to validate the model fit (convergence, Heywood cases, etc.).

---

## Repo layout

- `config/` : example config + (optional) model text file
- `bin/` : submission + helper scripts
- `slurm/` : sbatch scripts
- `scripts/` : R + bash pipeline logic
- `templates/multivariate_gwas_report/` : notebook template + expected folder skeleton
- `examples/` : non-pipeline examples (e.g., path diagrams)
