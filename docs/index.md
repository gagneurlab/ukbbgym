# Welcome to UKB-GYM

UKB-GYM provides a suite of DNAnexus applets and utility scripts for benchmarking variant annotation scoring models against LOFTEE correlations in quantitative traits and Olink proteomics. 

## Overview

Variant effect predictors are vital for interpreting genetic variants. UKB-GYM provides a standardized pipeline to compare these predictors using large-scale association data from the UK Biobank (UKB). 

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
