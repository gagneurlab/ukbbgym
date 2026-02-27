# Odds Ratio Computation and Bootstrapping Logic

## Overview

We ask: **are variants in annotated regions (e.g. frameshift, missense, splice) enriched among those with extreme phenotype effects?**

For each annotation type, we compute an odds ratio (OR) at increasingly extreme phenotype quantile thresholds, producing a curve that shows how enrichment grows with phenotype extremity.

---

## 1. Setup: Gene-Trait Pairs and Direction Correction

Each gene (region) is paired with its most strongly associated trait (670 gene-trait pairs, ~96 unique traits). The association direction (`loftee_corr_dir`) tells us whether loss-of-function variants increase or decrease the trait.

For each variant, we have `mean_pheno_value` — the average phenotype value across carriers of that variant. We assume this is drawn from a **standard Gaussian (z-score)**.

We apply direction correction so that "extreme" always means the same thing across gene-trait pairs:

```python
gp_base = (
    appv
    .join(id_region, on='id', how='inner')
    .join(gene_trait_df.lazy(), on=['region', 'phenotype'], how='inner')
    .with_columns(
        mean_pheno_value_dircor = pl.when(pl.col('loftee_corr_dir') == -1)
            .then(-pl.col('mean_pheno_value'))
            .otherwise(pl.col('mean_pheno_value'))
    )
    .select(['id', 'region', 'mean_pheno_value_dircor'])
)
# One row per (variant, gene-trait pair). ~22M rows × 3 cols.
# Columns: id (variant), region (gene), mean_pheno_value_dircor (direction-corrected z-score)
```

After this, a **large positive z-score** always means "phenotypically extreme in the expected direction."

---

## 2. Quantile Cutoffs in Z-Score Space

We define 101 quantile cutoffs on the percentile scale (log-spaced from 0 to 1 - 1e-5), then convert them to **z-score thresholds** using the inverse normal CDF:

```python
cutoffs_sorted = ...  # percentile cutoffs: 0.0, 0.11, ..., 0.99, ..., 0.99999
cutoffs_sorted = cutoffs_sorted.with_columns(
    zscore_cutoff = pl.Series(stats.norm.ppf(cutoffs_sorted['pheno_cutoff'].to_numpy()))
)
```

| Percentile | Z-score |
|-----------|---------|
| 0.50      | 0.00    |
| 0.90      | 1.28    |
| 0.95      | 1.64    |
| 0.99      | 2.33    |
| 0.999     | 3.09    |
| 0.99999   | 4.27    |

The **comparison** (is a variant "extreme"?) uses z-scores.
The **plot x-axis** shows the corresponding percentile quantiles (more interpretable).

---

## 3. The 2x2 Contingency Table

At each z-score threshold `z_c`, for each annotation, we classify every variant-gene-trait observation into a 2x2 table:

|                        | Annotated (e.g. frameshift) | Not annotated |
|------------------------|-----------------------------|---------------|
| **Extreme** (z > z_c)  | a                           | b             |
| **Not extreme** (z ≤ z_c) | c                        | d             |

The odds ratio is:

```
OR = (a / c) / (b / d) = (a * d) / (b * c)
```

- OR > 1 → annotated variants are enriched among extreme phenotype carriers
- OR = 1 → no enrichment

---

## 4. Efficient Computation via Cumulative Sums + Asof Join

Instead of cross-joining all observations with all 101 cutoffs, we use a cumsum approach:

### 4a. Group by z-score, count annotated vs not

For each (region, z-score value), count how many observations are annotated vs not:

```python
anno_data = (
    melted_anno
    .filter(pl.col('annotation') == annotation)
    .select(['id', 'region', 'annotation_score_dircor'])
)
# One row per (variant, gene). Columns: id, region, annotation_score_dircor (0 or 1).

region_zscore_counts = (
    gp_base.lazy()
    .join(anno_data.lazy(), on=['id', 'region'], how='inner')
    .group_by(['region', 'mean_pheno_value_dircor'])
    .agg(
        n_annotated = (pl.col('annotation_score_dircor') == 1).sum(),
        n_not_annotated = (pl.col('annotation_score_dircor') == 0).sum(),
    )
    .sort(['region', 'mean_pheno_value_dircor'])
    .with_columns(
        cum_annotated = pl.col('n_annotated').cum_sum().over('region'),
        cum_not_annotated = pl.col('n_not_annotated').cum_sum().over('region'),
    )
)
# One row per unique (region, z-score). Sorted by z-score within each region.
# Columns: region, mean_pheno_value_dircor, n_annotated, n_not_annotated,
#           cum_annotated (running sum ≤ this z-score), cum_not_annotated

region_totals = (
    region_zscore_counts
    .group_by('region')
    .agg(
        total_annotated = pl.col('n_annotated').sum(),
        total_not_annotated = pl.col('n_not_annotated').sum(),
    )
)
# One row per region. Columns: region, total_annotated, total_not_annotated.
```

### 4b. Asof join to map cutoffs to cumulative counts

For each (region, cutoff), find the cumulative count at the largest z-score ≤ the cutoff:

```python
per_region = (
    cutoffs_sorted
    .join(all_regions_df, how='cross')           # every (cutoff, region) pair
    .sort(['region', 'zscore_cutoff'])
    .join_asof(
        region_zscore_counts,
        left_on='zscore_cutoff',
        right_on='mean_pheno_value_dircor',
        by='region',
        strategy='backward',
    )
)
# One row per (region, cutoff). Shape: n_regions × n_cutoffs rows (e.g. 670 × 101 = 67,670).
# For each row, cum_annotated/cum_not_annotated = counts with z-score ≤ this cutoff.
```

### 4c. Derive the 2x2 cells

```python
a = n_ann_above   = total_annotated     - cum_annotated       # annotated & extreme
c = n_ann_below   = cum_annotated                             # annotated & not extreme
b = n_notann_above = total_not_annotated - cum_not_annotated  # not annotated & extreme
d = n_notann_below = cum_not_annotated                        # not annotated & not extreme
```

This gives us a 2x2 table for every **(region, annotation, cutoff)** combination.

---

## 5. Pre-Aggregation into a 4D Array

All per-region counts are stored in a numpy array for fast bootstrapping:

```python
counts_4d = np.zeros((n_annotations, n_regions, n_cutoffs, 4), dtype=np.int64)
# Shape: (7, 670, 101, 4) ≈ 15 MB
# Axis 0: annotation type
# Axis 1: gene-trait pair (region)
# Axis 2: quantile cutoff
# Axis 3: [a, c, b, d] = the four cells of the 2x2 contingency table
```

---

## 6. Point Estimate: Pooled OR

The point estimate sums counts across **all** gene-trait pairs into one pooled 2x2 table per (annotation, cutoff):

```python
total_counts = counts_4d.sum(axis=1)
# Shape: (n_annotations, n_cutoffs, 4). One pooled 2x2 table per (annotation, cutoff).

point_or = (total_counts[:, :, 0] / total_counts[:, :, 1]) / \
           (total_counts[:, :, 2] / total_counts[:, :, 3])
#         = (a_pooled / c_pooled) / (b_pooled / d_pooled)
# Shape: (n_annotations, n_cutoffs). One OR value per (annotation, cutoff).
```

This is a single OR computed from all 670 gene-trait pairs combined — **not** an average of per-pair ORs.

---

## 7. Bootstrap Confidence Intervals

### Why not Woolf CIs?

The standard Woolf SE formula assumes independent observations. But our ~670 gene-trait pairs share only ~96 unique traits (~7 genes per trait), so observations across genes with the same trait are correlated. Woolf CIs would be too narrow.

### Bootstrap procedure

We resample **gene-trait pairs** (regions) with replacement. Each bootstrap iteration:

1. Draw 670 regions with replacement
2. Sum their pre-computed counts → one resampled pooled 2x2 table
3. Compute the pooled OR from the resampled table

```python
n_boot = 1000
boot_ors = np.full((n_boot, n_annotations, n_cutoffs), np.nan)
# Shape: (1000, n_annotations, n_cutoffs). OR for each bootstrap iteration.

for b in range(n_boot):
    idx = rng.choice(n_regions, size=n_regions, replace=True)
    sampled = counts_4d[:, idx, :, :].sum(axis=1)
    # sampled shape: (n_annotations, n_cutoffs, 4). One resampled pooled 2x2 table.
    boot_ors[b] = (sampled[:, :, 0] / sampled[:, :, 1]) / \
                  (sampled[:, :, 2] / sampled[:, :, 3])
```

### Why this is fast

The expensive joins and cumsums happen once (Step 4). Each bootstrap iteration is just:
- **Fancy index**: `counts_4d[:, idx, :, :]` — select 670 rows (with repeats)
- **Sum**: `.sum(axis=1)` — collapse to pooled table
- **Divide**: compute OR

All annotations and cutoffs are processed simultaneously. 1000 iterations take seconds.

### Computing CIs

95% CIs use the **percentile method** directly on the bootstrap distribution of OR (no log transform needed, since percentiles are invariant to monotone transforms like exp/log):

```python
boot_ci_lower = np.nanpercentile(boot_ors, 2.5, axis=0)
boot_ci_upper = np.nanpercentile(boot_ors, 97.5, axis=0)
# Shape: (n_annotations, n_cutoffs). Lower/upper 95% CI bounds (in OR scale).
```

### Output DataFrame

```python
or_df = pl.DataFrame({
    'annotation': ...,            # annotation name (str)
    'pheno_cutoff': ...,          # percentile cutoff, for plotting (f64)
    'zscore_cutoff': ...,         # z-score threshold used for comparison (f64)
    'n_dis_above_cutoff': ...,    # a: annotated & extreme (i64)
    'n_notdis_above_cutoff': ..., # c: annotated & not extreme (i64)
    'n_dis_below_cutoff': ...,    # b: not annotated & extreme (i64)
    'n_notdis_below_cutoff': ..., # d: not annotated & not extreme (i64)
    'odds_ratio': ...,            # pooled OR = (a/c) / (b/d) (f64)
    'ci_lower': ...,              # bootstrap 2.5th percentile of OR (f64)
    'ci_upper': ...,              # bootstrap 97.5th percentile of OR (f64)
})
# One row per (annotation, cutoff). Filtered to finite OR > 0.
```

---

## 8. Visualization

The plot shows:
- **X-axis**: Extreme phenotype quantile (percentile scale, log-spaced) — e.g. "1e-2" = top 1%
- **Y-axis**: Odds ratio (log scale)
- **Lines**: One per annotation type
- **Ribbons**: 95% bootstrap CIs

The x-axis displays quantiles (not z-scores) because quantiles are more interpretable, even though the underlying comparison uses z-score thresholds.
