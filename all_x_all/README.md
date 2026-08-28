# All of Us All-by-All Summary Statistics for UKBBGym

## Scope

This replicates **only the correlation heatmap** of the [genebass](../genebass/README.md)
benchmark, against effect sizes from All of Us's [All by All](https://databrowser.researchallofus.org/genomic-variants)
(AbA) summary statistics instead of Genebass. It is not a general reimplementation — everything
about the method (APPV-as-rank-correlation-target, direction correction, the pairwise
significantly-beats heatmap and its topological-sort axis ordering) is the same as `../genebass/`;
see that README for the method itself. This directory documents only what differs: the data
source and what's implemented.

The one analysis notebook is [`analysis/correlations/pheno_correlations_META.ipynb`](analysis/correlations/pheno_correlations_META.ipynb),
which reads the AbA META (cross-ancestry meta-analysis) effect sizes and produces the pairwise
missense-predictor correlation heatmap plus its concordance with the full UKBBGym ranking.

## This does not run locally

Everything here — extraction and the analysis notebook alike — runs on the
[All of Us Researcher Workbench](https://www.researchallofus.org/), against data that cannot
leave it. `CONFIG_DIR` and `sys.path` in the analysis notebook resolve relative to wherever you
clone this repo *on the Workbench*; the extraction notebooks under `utils/` reference the
Workbench's own Hail/Spark cluster and Cloud Storage buckets and cannot be pointed elsewhere.

## `utils/` — AbA extraction (Hail/Spark, Workbench only)

Query and QC the All-by-All exome summary statistics before they reach the analysis notebook:

| Notebook | Does |
|---|---|
| [`query_rvat_exomes_AbA.ipynb`](utils/query_rvat_exomes_AbA.ipynb) | Rare-variant single-variant exome query |
| [`query_rvat_genes_AbA.ipynb`](utils/query_rvat_genes_AbA.ipynb) | Rare-variant gene-level (burden) query |
| [`query_exomes_META_AbA.ipynb`](utils/query_exomes_META_AbA.ipynb) | Cross-ancestry META exome query — this is what `pheno_correlations_META.ipynb` reads |
| [`extract_variant_metadata.ipynb`](utils/extract_variant_metadata.ipynb) | Variant-level metadata extraction |
| [`variant_qc_check.ipynb`](utils/variant_qc_check.ipynb) | QC checks on the extracted variant set |
| [`uncorrelated_phenos.ipynb`](utils/uncorrelated_phenos.ipynb) | Selects a near-independent phenotype subset |
| [`loftee_correlations.ipynb`](utils/loftee_correlations.ipynb) | LOFTEE-correlation-based association selection, matching `../ukbb/`'s association-filtering step (association selection, not a figure — hence `utils/` rather than `analysis/`) |

These are Hail/Spark notebooks and need compute sizing considerations similar to
[`genebass/utils/01_read_hail_sumstats.ipynb`](../genebass/utils/01_read_hail_sumstats.ipynb) —
run them on the Workbench's Spark cluster, not interactively on a small instance.

## `analysis/correlations/`

- [`pheno_correlations_META.ipynb`](analysis/correlations/pheno_correlations_META.ipynb) — the
  correlation heatmap. Reads [`../configs/all_x_all/`](../configs/all_x_all/) for the predictor
  set and variant-class filters (same YAML structure as `../configs/genebass/`, but scored
  against AbA's column names — e.g. `polyphen_score` rather than genebass's `polyphen`).
- [`dataframe_from_heatmap.py`](analysis/correlations/dataframe_from_heatmap.py) — helper the
  notebook imports for the pairwise heatmap's row/column tool ordering.

Figure output (`FIG_DIR`) resolves to `paper_figures/` at the repo root, same convention as
`../genebass/`.
