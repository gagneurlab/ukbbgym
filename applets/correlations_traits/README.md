# correlations_pheno

Computes Spearman correlations between variant annotation scores and phenotype means for benchmarking variant scoring models.

## Overview

This applet analyzes how well annotation scores (built-in tools like AlphaMissense, REVEL, SpliceAI, etc., plus optional custom models) correlate with average phenotype values for variants in a set of gene-trait associations.

Key features:
- Supports quantitative trait phenotypes (INT-transformed)
- Flexible annotation category selection (plof, missense, conservation, splicing, regulatory, etc.)
- Optional custom model benchmarking against built-in scores
- Variant class presets (missense, missense_structured, intron, enhancer, promoter, etc.) or custom filters
- Direction-corrected correlation accounting for LOFTEE effect direction
- Outputs: correlation results, boxplot, and pairwise Wilcoxon heatmap

## Inputs

| Input | Required | Type | Description |
|-------|----------|------|-------------|
| `associations_parquet` | Yes | file | Gene-trait associations with LOFTEE correlation. Columns: `region`, `phenotype`, `correlation`. |
| `annotations_parquet` | Yes | file | Variant annotations with all scores. Columns: `id`, `region`, + score columns. |
| `appv_parquet` | Yes | file | Average phenotype per variant (INT-transformed). Columns: `id`, `phenotype`, `mean_pheno_value`, `n_individuals`. |
| `custom_scores_parquet` | No | file | Custom model scores. Format: `id` + optional `region` + float columns per model. |
| `model_labels` | No | string | Comma-separated display names for custom score columns. |
| `model_directions` | No | string | Comma-separated `1` or `-1` per custom model (1=higher is worse, -1=lower is worse). Default: `1` for all. |
| `fill_missing_scores` | No | bool | If True, fill custom score nulls with 0.0. Default: False. |
| `annotation_categories` | No | string | Comma-separated categories to include. Default: `plof,missense,conservation,splicing,regulatory`. |
| `variant_class` | No | string | Preset name or `"other"`. Default: `all_variants`. |
| `variant_class_config_yaml` | No | file | Required if `variant_class="other"`. YAML with variant filtering rules. |
| `annotation_config_yaml` | No | file | Override bundled annotation config (label/color/direction per tool). |
| `mac_threshold` | No | int | Max `n_individuals` per variant. Default: 20. |
| `only_snps` | No | bool | Exclude indels. Default: True. |
| `only_clinvar` | No | bool | Include only ClinVar-annotated variants. Default: False. |
| `exclude_clinvar` | No | bool | Exclude ClinVar-annotated variants. Default: False. |

## Output

| Output | Name | Description |
|--------|------|-------------|
| `results_parquet` | file | Correlation table. Columns: `region`, `phenotype`, `annotation`, `n_variants`, `correlation`, `category`, `label`, `color`, `annotation_dir`, `corr_beta`, `corr_beta_rescaled`. |
| `boxplot_svg` | file | Boxplot: `corr_beta` per annotation, ordered by median. |
| `heatmap_svg` | file | Pairwise Wilcoxon heatmap: mean_diff per annotation pair. |

## Build

```bash
dx build applets/correlations_pheno/ --destination <project>:/applets/ -f
```

## Details

- Phenotypes are INT-transformed (z-scored) using Inverse Normal Transformation
- Spearman correlation computed via ranking with "average" tie handling
- Direction correction: `corr_beta = correlation × loftee_corr_dir × annotation_dir`
- Filtering: `n_variants > 100` per gene-trait-annotation triple
- Instance: `mem3_ssd1_v2_x16` (128 GB). Runtime ~40 minutes.
