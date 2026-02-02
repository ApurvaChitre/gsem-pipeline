# Environment

If you already ran the MPH pipeline, you likely already have the conda environment (e.g. `genomic_sem`) that contains:

- R + required CRAN/Bioconductor packages
- GenomicSEM
- PLINK

This repo includes `env/environment.example.yml` as a *starting point*, but for maximum reproducibility you should use the **lock file** from the MPH pipeline repo (if available).

If `library(GenomicSEM)` fails inside your conda env, you can try:

```bash
Rscript env/install_genomicsem_from_github.R
```
