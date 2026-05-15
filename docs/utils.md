# Utilities & Scripts

The `utils/` directory contains a variety of scripts and Jupyter Notebooks used internally to process associations, format BCFs to Parquet, and prepare data.

## 1_REGENIE_rvat

This directory contains resources for running REGENIE for rare variant association testing. It includes configuration templates and submission scripts required to properly setup and execute distributed jobs across the UK Biobank dataset.

## 2_association_filtering

Once associations are calculated (e.g., via REGENIE), they must be filtered down to high-confidence sets before being used in downstream analysis.
The notebooks in this folder handle:
- Filtering associations based on p-value thresholds.
- Extracting effect directions (e.g., LOFTEE effect directions).
- Consolidating the final gene-trait pairs for the `associations_parquet` inputs.

## bcf2parquet & annotations

Before variant scores can be analyzed, the raw UKB genetic data (often stored as BCF files) must be joined with consequence and scoring annotations.
- The `bcf2parquet` scripts handle the conversion of binary format files into the optimized Parquet structures required by the applets.
- The `annotations` notebooks add standard UKB-GYM variant annotations (e.g., from VEP, structural prediction models, or conservation scores) to the resulting Parquet files.
