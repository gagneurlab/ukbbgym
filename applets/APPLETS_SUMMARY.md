# DNAnexus Analysis Applets for Variant Annotation Benchmarking

## Overview

Six published DNAnexus applets for benchmarking variant annotation scoring models against LOFTEE correlations in quantitative traits and Olink proteomics, as well as for preparing variant-level phenotype aggregates.

## Applets

### 1. `avg_pheno_per_variant_traits`
**Average Phenotype Per Variant for quantitative traits**

Computes the aggregated mean phenotype value and individual count per variant from the initial processed UK Biobank parquet files.

**Key Inputs:** Initial association parquets from upstream processing
**Key Outputs:** `appv_parquet` containing `id`, `phenotype`, `mean_pheno_value`, and `n_individuals`.

---

### 2. `avg_pheno_per_variant_olink`
**Average Phenotype Per Variant for Olink proteomics**

Identical logic to `avg_pheno_per_variant_traits` but designed for PROTRIDER-corrected protein abundance data.

---

### 3. `correlations_traits`
**Spearman correlations for quantitative trait phenotypes**

Computes how well each annotation score correlates with average phenotype values for variants in a gene-trait association set.

**Workflow:**
1. Load pre-filtered gene-trait associations with LOFTEE effect directions
2. Filter and select variant annotations by class and category
3. Merge optional custom model scores
4. Melt annotations to long format and join with APPV data
5. Compute Spearman correlations (via ranking) per gene-trait-annotation
6. Direction-correct: `corr_beta = correlation × loftee_corr_dir × annotation_dir`
7. Output results, boxplot, and pairwise Wilcoxon heatmap

**Key Inputs:** associations, annotations, APPV parquets; optional custom scores
**Key Parameters:** `annotation_categories`, `variant_class`, `mac_threshold`, `only_snps`, ClinVar filters
**Key Outputs:** correlation results parquet, boxplot SVG, heatmap SVG

---

### 4. `correlations_olink`
**Spearman correlations for Olink proteomics**

Identical logic to `correlations_traits` but for PROTRIDER-corrected protein abundance data (no INT transformation).

**Differences from traits:**
- APPV phenotypes are already normalized (PROTRIDER), not INT-transformed
- All else identical

---

### 5. `top_n_vars_traits`
**Top-N variant benchmarking for quantitative traits**

Ranks variants by annotation scores within each gene and computes cumulative phenotype effects at each rank, then compares across annotations.

**Workflow:**
1. Load associations, annotations, APPV (same setup as correlations)
2. Direction-correct phenotype values by LOFTEE effect
3. Rank variants within (annotation, gene) by score descending
4. For each rank 1 to max_gene_rank: compute mean phenotype at that rank within each gene
5. Cumulative mean: `cum_mean = cum_sum / ordinal_rank` per gene
6. Aggregate across genes: mean, std, SEM, CI for each (annotation, rank)
7. Compute per-gene top-N means for each window size in `target_k`
8. Pairwise Wilcoxon signed-rank tests for annotation pairs, faceted by top-N window
9. Output cumulative stats, per-gene top-N, line plot, and heatmap

**Key Inputs:** associations, annotations, APPV parquets; optional custom scores
**Key Parameters:** `target_k` (comma-separated top-N values, default "5,10"), `max_gene_rank` (default 50)
**Key Outputs:** cumulative stats parquet, per-gene top-N parquet, line plot SVG, faceted heatmap SVG

---

### 6. `top_n_vars_olink`
**Top-N variant benchmarking for Olink proteomics**

Identical logic to `top_n_vars_traits` but for PROTRIDER-corrected protein data.

---

## Shared Input Specification (Analysis Applets)

The analysis applets (`correlations` and `top_n_vars`) accept the same core inputs:

| Input | Required | Type | Description |
|-------|----------|------|-------------|
| `associations_parquet` | Yes | file | Gene-trait/protein associations with LOFTEE correlation direction. Columns: `region`, `phenotype`, `correlation`. |
| `annotations_parquet` | Yes | file | Variant annotations with all scoring tools. Columns: `id`, `region`, + score columns. |
| `appv_parquet` | Yes | file | Average phenotype/protein per variant. Columns: `id`, `phenotype`, `mean_pheno_value`, `n_individuals`. |
| `custom_scores_parquet` | No | file | Custom model scores. Format: `id` + optional `region` + one float column per model. |
| `model_labels` | No | string | Comma-separated display names for custom model columns, in column order. E.g., `"MyModel v1,MyModel v2"`. |
| `model_directions` | No | string | Comma-separated `1` or `-1` per custom model, in column order. `1` = higher is worse, `-1` = lower is worse. Default: `1` for all. |
| `fill_missing_scores` | No | bool | If True, fill custom score nulls with 0.0 (neutral). If False, exclude variants without a custom score. Default: False. |
| `annotation_categories` | No | string | Comma-separated annotation categories to include. Options: `plof`, `missense`, `conservation`, `splicing`, `regulatory`, `genetic_diversity`, `plof_consequences`, `clinvar`, `protein_domains`, `5utr`, `3utr`, `indel_lengths`, `non_coding_regions`, etc. Default: `plof,missense,conservation,splicing,regulatory`. |
| `variant_class` | No | string | Preset variant class name or `"other"`. Options: `all_variants`, `all_coding`, `all_non_coding`, `missense`, `missense_structured`, `missense_non_structured`, `missense_disordered`, `missense_lip`, `intron`, `enhancer_encode`, `promoter_encode`, `tf_encode`, `proximal_promoter`, `core_promoter`. Default: `all_variants`. |
| `variant_class_config_yaml` | No | file | Required if `variant_class="other"`. User-provided YAML with same structure as bundled `config_variant_classes.yaml`. |
| `annotation_config_yaml` | No | file | Override bundled annotation config. YAML with `rare_variant_annotations` dict mapping category → annotation → {`color`, `label`, optional `direction`}. |
| `mac_threshold` | No | int | Max `n_individuals` per variant. Default: 20. |
| `only_snps` | No | bool | Exclude indels (ref/alt length > 1). Default: True. |
| `only_clinvar` | No | bool | Include only ClinVar-annotated variants. Default: False. |
| `exclude_clinvar` | No | bool | Exclude ClinVar-annotated variants. Default: False. |

### top_n_vars-only additional inputs

| Input | Type | Description |
|-------|------|-------------|
| `target_k` | string | Comma-separated top-N window sizes for heatmap faceting. Default: `"5,10"`. |
| `max_gene_rank` | int | Max rank depth for cumulative z-score line plot. Default: 50. |

## Custom Scores Format

Users can upload custom model scores as a single parquet file with:
- **Required column:** `id` (variant identifier, matching annotations/APPV)
- **Optional column:** `region` (gene/region identifier)
  - If present, custom scores are joined on `[id, region]`
  - If absent, joined on `[id]` only
- **Score columns:** All remaining float columns are treated as model scores
  - Column names become the model identifiers
  - Supply `model_labels` and `model_directions` to customize display

**Example:**
```
id                   region          mymodel_v1  mymodel_v2
chr1:100:A:T         ENSG0000001     0.75        0.82
chr1:101:C:G         ENSG0000001     0.45        0.51
...
```

## Configuration Files

### `config_correlations.yaml` (for correlations applets)
Defines annotation metadata for correlation-based analysis. Used by correlations_traits/olink.
- Keys: annotation categories (plof, missense, conservation, splicing, regulatory, etc.)
- Per-annotation: `color` (hex), `label` (display name), optional `direction` (1 or -1)

### `config_odds.yaml` (for top-N applets)
Defines annotation metadata for top-N analysis. Used by top_n_vars_traits/olink.
- Same structure as config_correlations.yaml but with different color/label choices

### `config_variant_classes.yaml` (shared)
Defines preset variant filtering classes. Used by all analysis applets.
- Keys: class names (missense, intron, enhancer_encode, etc.)
- Per-class: `variant_filtering` (list of Polars filter expressions as strings), `tool_categories` (which to include), `x_label` (display name)

All configs are bundled in each applet; users can override via optional YAML inputs.

## Building and Running

### Build an applet
Note: All applets are publicly available. However, if building from source:
```bash
dx build applets/correlations_traits/ --destination project-XXX:/applets/ -f
```

### Run via CLI
```bash
dx run correlations_traits \
  -i associations_parquet=file-XXX \
  -i annotations_parquet=file-YYY \
  -i appv_parquet=file-ZZZ \
  -i variant_class="missense" \
  -i annotation_categories="missense,conservation" \
  --destination project-XXX:/outputs/ \
  --brief
```

### Run with custom scores
```bash
dx run top_n_vars_traits \
  -i associations_parquet=file-XXX \
  -i annotations_parquet=file-YYY \
  -i appv_parquet=file-ZZZ \
  -i custom_scores_parquet=file-CUSTOM \
  -i model_labels="MyModel v1,MyModel v2" \
  -i model_directions="1,-1" \
  -i target_k="5,10,20" \
  -i max_gene_rank=100 \
  --destination project-XXX:/outputs/ \
  --brief
```

## Outputs

### Correlations applets
- **results_parquet:** Full correlation table with all computed metrics
- **boxplot_svg:** Boxplot of correlation_beta per annotation, ordered by median
- **heatmap_svg:** Pairwise Wilcoxon heatmap (mean_diff with significance stars)

### Top-N applets
- **cum_stats_parquet:** Cumulative z-score statistics by rank (mean, std, SEM, 95% CI)
- **gene_topn_parquet:** Per-gene top-N z-scores for each target_k window
- **lineplot_svg:** Cumulative mean z-score by rank with ±1 SEM ribbon per annotation
- **heatmap_svg:** Pairwise Wilcoxon heatmap faceted by target_k window (mean_diff with significance stars)

## Instance Requirements

All applets run on `mem3_ssd1_v2_x16` (128 GB) in eu-west-2.

**Approximate runtime:**
- correlations_traits/olink: ~40 minutes
- top_n_vars_traits/olink: ~30 minutes
- avg_pheno_per_variant: ~1-3 hours (depending on data scale)

(Runtime depends on number of genes, variants, and annotations selected)

## Key Design Decisions

1. **No LOEUF stratification:** Per-gene LOEUF quartile plots are not included in applet outputs (available in Jupyter notebooks)
2. **SVG outputs only:** Plots are exported as SVG for publication quality and easy inclusion in reports
3. **Flexible annotation selection:** Users can mix built-in and custom tools in a single run
4. **Pre-filtered associations:** Applets assume associations are already filtered by p-value; they use all provided gene-trait pairs
5. **Custom scoring metadata:** Users specify model names, directions, and missing-value handling at run time (not in config files)
6. **Publishable design:** No hard-coded project paths or local file references; all inputs are explicit file uploads

## Notes for Users

- **Associations file:** Must contain merged gene-trait pairs and LOFTEE correlations in columns `region`, `phenotype`, `correlation`
- **Custom score nulls:** By default, variants without a custom score are excluded for that model. Set `fill_missing_scores=True` to impute nulls with 0.0
- **Direction conventions:**
  - `annotation_dir=1` (default): higher score = more deleterious
  - `annotation_dir=-1`: lower score = more deleterious
  - Custom models default to `direction=1`; override via `model_directions` input
- **MAC threshold:** Variants with `n_individuals > mac_threshold` are excluded from correlations (default 20)
- **Variant filtering:** By default, only SNPs are included; set `only_snps=False` to include indels
