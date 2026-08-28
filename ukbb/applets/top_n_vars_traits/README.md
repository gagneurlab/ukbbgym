# top_n_vars_pheno

Benchmarks variant scoring models by comparing cumulative phenotype z-scores when ranking variants by model scores (top-N analysis).

## Overview

This applet ranks variants by annotation scores within each gene and computes cumulative phenotype effects. It evaluates whether higher-scoring variants (according to each annotation) tend to have stronger phenotypic effects, comparing built-in scores against optional custom models.

Key features:
- Per-gene top-N rankings: variants ranked by score, phenotypes averaged over top N
- Cumulative z-score statistics: per-gene cumulative mean at each rank (1 to max_gene_rank)
- Pairwise comparisons: Wilcoxon signed-rank tests for each annotation pair
- Multiple window sizes: compare performance at different top-N thresholds
- Supports flexible annotation categories and variant class presets
- Outputs: cumulative stats, per-gene top-N means, line plot, and faceted heatmap

## Inputs

| Input | Required | Type | Description |
|-------|----------|------|-------------|
| `associations_parquet` | Yes | file | Gene-trait associations with LOFTEE correlation. Columns: `region`, `phenotype`, `correlation`. |
| `annotations_parquet` | Yes | file | Variant annotations with all scores. Columns: `id`, `region`, + score columns. |
| `appv_parquet` | Yes | file | Average phenotype per variant (INT-transformed). Columns: `id`, `phenotype`, `mean_pheno_value`, `n_individuals`. |
| `custom_scores_parquet` | No | file | Custom model scores. Format: `id` + optional `region` + float columns per model. |
| `model_labels` | No | string | Comma-separated display names for custom score columns. |
| `model_directions` | No | string | Comma-separated `1` or `-1` per custom model. Default: `1` for all. |
| `fill_missing_scores` | No | bool | If True, fill custom score nulls with 0.0. Default: False. |
| `annotation_categories` | No | string | Comma-separated categories. Default: `plof,missense,conservation,splicing,regulatory`. |
| `variant_class` | No | string | Preset name or `"other"`. Default: `all_variants`. |
| `variant_class_config_yaml` | No | file | Required if `variant_class="other"`. |
| `annotation_config_yaml` | No | file | Override bundled annotation config. |
| `mac_threshold` | No | int | Max `n_individuals` per variant. Default: 20. |
| `only_snps` | No | bool | Exclude indels. Default: True. |
| `only_clinvar` | No | bool | Include only ClinVar-annotated variants. Default: False. |
| `exclude_clinvar` | No | bool | Exclude ClinVar-annotated variants. Default: False. |
| `target_k` | No | string | Comma-separated top-N window sizes. Default: `5,10`. |
| `max_gene_rank` | No | int | Max rank depth for line plot. Default: 50. |

## Output

| Output | Name | Description |
|--------|------|-------------|
| `cum_stats_parquet` | file | Cumulative z-score by rank. Columns: `annotation`, `ordinal_rank`, `mean_cum_zscore`, `std_cum_zscore`, `n_genes`, `se_cum_zscore`, `ci_lower`, `ci_upper`, `label`, `color`. |
| `gene_topn_parquet` | file | Per-gene top-N means. Columns: `region`, `annotation`, `gene_topn_zscore`, `n_variants_used`, `top_n`. |
| `lineplot_svg` | file | Cumulative z-score line plot with ±1 SEM ribbons. |
| `heatmap_svg` | file | Pairwise Wilcoxon heatmap, faceted by top-N window. |

## Build

```bash
dx build applets/top_n_vars_pheno/ --destination <project>:/applets/ -f
```

## Details

- Within each (annotation, gene), variants are ranked by score descending (method='max' for ties)
- Per-rank z-scores are averaged over tied variants (ordinal rank fills gaps)
- Cumulative mean: `cum_sum(z) / ordinal_rank` within each (annotation, gene)
- Pairwise comparisons: paired Wilcoxon test across genes for each annotation pair
- Instance: `mem3_ssd1_v2_x16` (128 GB). Runtime ~30 minutes.
