# Welcome to UKBBGym

UK Biobank Gym (UKBBGym) is a reusable population-scale benchmarking framework that evaluates variant deleteriousness predictors against human population phenotypes. By anchoring evaluation in observed per-variant phenotypic effects, spanning 670 gene-trait associations and 1,076 gene-protein abundance associations across over 40 million naturally occurring rare variants in the UK Biobank, UKBBGym avoids the ascertainment biases of clinical-label benchmarks. It provides a standardized set of associations, metrics, and reproducible workflows on the UK Biobank Research Analysis Platform (RAP) to directly compare computational scoring methods and experimental assays on a shared, population-derived ground truth.

```{note} UKBBGym Summary Statistics
UKBBGym's benchmarks are also replicated against externally available summary statistics cohorts (All of Us AllxAll, Genebass), for settings where only gene/variant-level summary statistics are available rather than individual-level UK Biobank data. See the separate **UKBBGym Summary Statistics documentation** (link to be added once hosting is finalized; `ukbgym_sumstats` is currently a private repository).
```

## Overview

Variant effect predictors are vital for interpreting genetic variants. UKBBGym provides a standardized pipeline to compare these predictors using large-scale association data from the UK Biobank (UKB). To guard against confounding associations, variants evaluated in UKBBGym do not occur in more than 20 carriers.

The workflow is divided into two major phases:
1. **Data Processing:** Processing UKB raw data into structured parquets and aggregating phenotype effects per variant (`avg_pheno_per_variant` applets).
2. **Analysis:** Computing correlations and top-N variant statistics across different functional annotations to benchmark scoring models (`correlations` and `top_n_vars` applets).

## Prerequisites

To use the tools provided in this repository, you will need:
- Access to the [UK Biobank Research Analysis Platform (UKB-RAP)](https://ukbiobank.dnanexus.com/).
- Basic understanding of Python and working with Jupyter Notebooks.
- `dx-toolkit` installed locally (if running commands via CLI).

## Getting Started

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ShubhankarLondhe/ukbgym.git
   cd ukbgym
   ```

2. **Setup your environment:**
   We highly recommend using `pre-commit` to maintain code quality.
   ```bash
   pip install pre-commit
   pre-commit install
   ```

3. **Accessing the Applets:**
   The DNAnexus applets used for processing and analysis are already compiled and publicly available on the UKB-RAP. You do **not** need to build them from source to use them. See the subsequent sections for details on running these applets via the GUI or CLI.

```{toctree}
:hidden:
:maxdepth: 2

data_processing
analysis
utils
config
```
