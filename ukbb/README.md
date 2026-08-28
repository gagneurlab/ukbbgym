# UKBBGym — UK Biobank RAP pipeline

This is the **primary UKBBGym benchmark**: it runs on individual-level UK Biobank data on the
UK Biobank Research Analysis Platform (RAP). It is the source of the manuscript's main figures
(1, 2, 3, 5a–5d, 6a and most supplementary figures) and of the association files the two
summary-statistics reconstructions (`../genebass/`, `../all_x_all/`) are built against.

**None of this can be run or tested outside the RAP** — individual-level UK Biobank data cannot
leave it. This directory is a reference for RAP users and a record of how the published results
were produced, not a locally runnable pipeline. If you don't have RAP access, see
`../genebass/README.md` for the setup that runs on a laptop.

## Overview

UKBBGym provides a standardized pipeline to compare variant effect predictors using large-scale
association data from the UK Biobank. To guard against confounding associations, variants
evaluated do not occur in more than 20 carriers (allele count ≤ 20).

The workflow has two phases:
1. **Data processing** — raw UKB data → structured parquets → average phenotype per variant
   (`avg_pheno_per_variant` applets).
2. **Analysis** — correlations and top-N variant statistics across functional annotations
   (`correlations` and `top_n_vars` applets).

## Prerequisites

- Access to the [UK Biobank Research Analysis Platform (UKB-RAP)](https://ukbiobank.dnanexus.com/).
- Basic familiarity with Python and Jupyter notebooks.
- `dx-toolkit` installed locally if you want to run applets from the CLI rather than the GUI.

## 1. Data processing

Raw UK Biobank BGEN/VCF and phenotype data are first converted into optimized Parquet files by
an upstream pipeline (maintained outside this repository); the output is a set of association
and basic-annotation parquets. From there:

- [`utils/bcf2parquet/`](utils/bcf2parquet/) — BCF → Parquet conversion notebooks (variant
  file-list construction, parquet consolidation, variant metadata, EUR subsetting).
- [`utils/3_average_pheno_per_variant/`](utils/3_average_pheno_per_variant/) — builds the
  average-phenotype-per-variant (APPV) tables locally, matching what the RAP applets below
  compute.

On the RAP itself, APPV is computed by two published applets:
- `avg_pheno_per_variant_traits` — quantitative trait phenotypes.
- `avg_pheno_per_variant_olink` — Olink proteomics (PROTRIDER-corrected).

Both output an `appv_parquet` with `id`, `phenotype`, `mean_pheno_value`, `n_individuals`.

**Running via the RAP GUI:** open the applet in your UKB-RAP project → **Run** → attach the
input parquet(s) → set parameters (e.g. minimum individuals per variant) → launch.

**Running via the CLI:**
```bash
dx run avg_pheno_per_variant_traits \
  -i input_parquet=project-XXX:/path/to/associations.parquet \
  -i min_individuals=20 \
  --destination project-XXX:/outputs/appv_traits/ \
  --brief
```

Gene–trait associations themselves are selected upstream by REGENIE burden testing and
filtering — see [`../utils/REGENIE_burden_test/`](../utils/REGENIE_burden_test/) for the
association-testing scripts and the p-value/effect-direction filtering notebooks that produce
the `associations_parquet` these applets read.

## 2. Analysis applets

Four published DNAnexus applets benchmark variant annotation scores against the APPV tables —
see [`applets/APPLETS_SUMMARY.md`](applets/APPLETS_SUMMARY.md) for full per-applet detail
(inputs, outputs, workflow steps). In short:

| Applet | Computes |
|---|---|
| `correlations_traits` | Spearman correlation between annotation score and APPV, quantitative traits |
| `correlations_olink` | same, Olink proteomics |
| `top_n_vars_traits` | rank variants by score, cumulative phenotype effect per rank, quantitative traits |
| `top_n_vars_olink` | same, Olink proteomics |

**Required inputs** (all four applets): `associations_parquet` (gene–trait/protein associations
with LOFTEE effect directions), `annotations_parquet` (variant annotation scores), `appv_parquet`
(from step 1 above).

**Running via the RAP GUI:** open the applet → **Run** → attach the three parquet inputs → set
options (`variant_class`, `annotation_categories`, `mac_threshold`, …) → launch.

**CLI example (`correlations_traits`):**
```bash
dx run correlations_traits \
  -i associations_parquet=project-XXX:/path/associations.parquet \
  -i annotations_parquet=project-XXX:/path/annotations.parquet \
  -i appv_parquet=project-XXX:/path/appv_traits.parquet \
  -i variant_class="missense" \
  -i annotation_categories="missense,conservation" \
  --destination project-XXX:/outputs/correlations/ \
  --brief
```

### Adding custom scores

Benchmark your own scoring model against the established tools by uploading it as a single
Parquet file: an `id` column matching the variant IDs in your annotations/APPV files, plus any
number of float columns treated as model scores.

```bash
dx run top_n_vars_traits \
  -i associations_parquet=project-XXX:/path/associations.parquet \
  -i annotations_parquet=project-XXX:/path/annotations.parquet \
  -i appv_parquet=project-XXX:/path/appv_traits.parquet \
  -i custom_scores_parquet=project-XXX:/path/custom_scores.parquet \
  -i model_labels="MyCustomModel v1,MyCustomModel v2" \
  -i model_directions="1,-1" \
  --destination project-XXX:/outputs/top_n/
```

### Outputs

`results_parquet` / `cum_stats_parquet` (all computed metrics, z-scores, standard errors) plus
publish-ready SVG plots — boxplots, line plots, pairwise Wilcoxon heatmaps.

## Configuration

Applets are configured through YAML files in [`../configs/ukbb/`](../configs/ukbb/); each applet
ships with defaults you can override via CLI/GUI input.

- **`config_variant_classes.yaml`** — the variant-filtering presets available via
  `-i variant_class=`. Each entry: `variant_filtering` (Polars filter expressions), 
  `tool_categories`, `x_label`. Provide your own class by passing
  `-i variant_class="other"` with a custom `variant_class_config_yaml`.
- **`config_correlations.yaml`** — visual metadata for `correlations_traits`/`correlations_olink`:
  per-tool `color`, `label`, and optional `direction` (±1, standardizing which way "more
  deleterious" points).
- **`config_odds.yaml`** — same structure, for `top_n_vars_traits`/`top_n_vars_olink` (distinct
  palette/labels tuned for cumulative z-score plots and heatmaps).
- **`config_odds_categories.yaml`**, **`config_expAssays.yaml`** — category-color mapping and
  the experimental-assay (SGE/MaveDB) predictor set, respectively.

## Local analysis notebooks

[`analysis/`](analysis/) holds the notebooks used to produce the manuscript's UKB-derived
figures directly from applet outputs (average-z-score facets, correlations, top-N genewise
statistics, ClinVar/LOEUF/protein-domain benchmarks) — these read applet output parquets and are
not runnable without them.

## Known gaps

The applet `README.md`/`dxapp.json` files under `applets/`, and a handful of notebooks under
`analysis/other_benchmarks/` (e.g. `clinvar_auroc_heatmap.ipynb`, `num_variants_in_other_benchmarks.ipynb`)
still hardcode paths from the pre-restructure repository (`CFG`, `FIG_DIR` pointing at the old
cluster layout), and four `analysis/` notebooks (`scatter_plot_1gene_olink.ipynb`,
`scatter_plot_1gene_pheno.ipynb`, `avg_zscore/avg_pheno_zscore_binary_annos.ipynb`,
`avg_zscore/avg_olink_zscore_binary_annos.ipynb`) reference a `config_pheno.yaml` that was
renamed away years ago and no longer exists. Since none of this runs outside the RAP/cluster
anyway, none of it was rewired in this pass — flagged here as a known follow-up.
