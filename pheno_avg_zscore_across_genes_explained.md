# How `pheno_avg_zscore_across_genes.ipynb` Computes Average Z-Scores in Sliding Windows

## Overview

This notebook computes the **mean direction-corrected phenotype z-score** across rare variants, ranked by various missense pathogenicity annotations, using a sliding window approach with bootstrap confidence intervals. The goal is to show that variants ranked as more deleterious by a given annotation tool (e.g. AlphaMissense, REVEL) also tend to have larger phenotypic effects.

---

## 1. Joins: How the Data Are Linked

The pipeline connects three core datasets through a series of joins:

### 1a. Gene-Trait Associations (`gene_trait_df`)

Significant gene-trait associations (FDR <= 0.05) are loaded and joined with LOFTEE correlation estimates. Each gene is kept only once — the association with the **largest absolute correlation** — to yield a single "most associated phenotype" per gene along with its correlation direction.

```python
gene_trait_df = (
    pl.read_parquet(f'{LOCAL_DIR}/{ASSOC_FILE}')
    .filter(pl.col('pval_fdr')<=0.05)
    .select(['region', 'phenotype', 'pval_fdr'])
)

# Join with LOFTEE correlations
gene_trait_df = (
    gene_trait_df
    .join(loftee_corrs, on=['region', 'phenotype'], how='inner')
    .drop_nans()
    .sort('loftee_corr_abs', descending=True)
    .unique(subset=["region"], keep="first", maintain_order=True)
)
```

This produces 670 unique gene-trait pairs, each with a `loftee_corr_dir` (+1 or -1) indicating whether LoF variants increase or decrease the phenotype.

### 1b. Variant Annotations (`melted_anno`)

The annotation file is filtered to the relevant gene regions, then restricted to **missense SNPs in LIP (linear interacting peptide) domains**. Selected annotation columns (e.g. `am_pathogenicity`, `revel_score`) are unpivoted from wide to long format, joined with the annotation config to get each annotation's expected direction, and direction-corrected so that higher rank always means "more deleterious". Variants with imputed/fill-NA scores are removed via an `anti` join.

```python
melted_anno = (
    anno.lazy()
    .unpivot(
        index=["id", "region"],
        on=selected_annos,
        variable_name="annotation",
        value_name="annotation_score"
    )
    # Join annotation config to get direction
    .join(
        anno_config_df.select(["annotation", "category", "annotation_dir"]).lazy(),
        on="annotation",
        how="left"
    )
    # Remove variants whose score was imputed (fill-NA)
    .join(
        anno_fillna_melted,
        on=['id', 'region', 'annotation'],
        how='anti'
    )
    # Direction-correct and rank
    .with_columns(
        annotation_score_dircor = pl.col('annotation_score') * pl.col("annotation_dir").cast(pl.Float32)
    )
    .with_columns(
        annotation_score_dircor_rank_desc = pl.col('annotation_score_dircor')
            .rank(method="max", descending=True).over(["annotation"]).cast(pl.Float32)
    )
    .collect(engine='streaming')
)
```

**Key detail on the anti-join**: Each annotation has a companion `<annotation>_is_na` column in the source parquet. The notebook unpivots these flags, filters to rows where the flag is 1 (i.e. the score was imputed), and uses an `anti` join on `(id, region, annotation)` to exclude those variant-annotation pairs from downstream analysis.

### 1c. Phenotype Average Per Variant (`pheno_appv`)

Pre-computed per-variant average phenotype values (mean phenotype z-score among carriers) are loaded, restricted to variants in the annotation set via a `semi` join, and further filtered to low MAC (minor allele count <= 20).

```python
anno_keys = anno.select(pl.col('id').unique()).lazy()

pheno_appv = (
    pheno_appv
    .join(anno_keys, on='id', how='semi')
    .filter(pl.col('n_individuals') <= mac)
    .select(['id', 'phenotype', 'mean_pheno_value', 'n_individuals'])
)
```

### 1d. The Final Join: Variant Z-Scores + Annotation Ranks

The three datasets are combined:

1. `pheno_appv` is joined with `id_region` (from the annotation table) on `id` to attach the gene (`region`).
2. Then joined with `gene_trait_df` on `(region, phenotype)` to get only the most significant gene-trait phenotype and its LOFTEE correlation direction.
3. The `mean_pheno_value` is **direction-corrected** by multiplying by `loftee_corr_dir`, so that deleterious effects always appear as positive z-scores.
4. Finally, this is joined with the melted annotation ranks on `(id, region)` to pair each variant's z-score with its rank under each annotation.

```python
variant_zscores = (
    pheno_appv
    .join(id_region, on='id', how='inner')
    .join(gene_trait_df.lazy(), on=['region', 'phenotype'], how='inner')
    .with_columns(
        mean_pheno_value = pl.col('mean_pheno_value') * pl.col('loftee_corr_dir').cast(pl.Float32)
    )
    .select(['id', 'region', 'mean_pheno_value'])
    .collect(engine='streaming')
)

ranked_zscores = (
    variant_zscores.lazy()
    .join(
        melted_anno.lazy().select(['id', 'region', 'annotation', 'annotation_score_dircor_rank_desc']),
        on=['id', 'region'],
        how='inner'
    )
    .sort(['annotation', 'annotation_score_dircor_rank_desc'])
    .with_columns(
        row_pos = pl.col('annotation').cum_count().over('annotation'),
    )
    .collect(engine='streaming')
)
```

---

## 2. Sliding Windows: How They Are Computed

### Parameters

```python
window_size = 1_000   # Number of variants in each window
step_size = 10        # Step between consecutive window positions
```

### Bin Endpoints

For each annotation, variants are sorted by their direction-corrected rank (most deleterious first). The **bin endpoints** are the right edges of each sliding window, generated as:

```python
bin_ends = np.arange(window_size, N + 1, step_size)
```

So the first window covers variants at positions `[0, 1000)`, the next covers `[10, 1010)`, then `[20, 1020)`, etc. Each window always contains exactly `window_size` variants.

### Efficient Computation via Prefix Sums

Rather than recomputing the mean for each window from scratch, the notebook uses **prefix sums** for O(1) per-window computation:

```python
def sliding_window_means(zscores, bin_ends, window_size, weights=None):
    if weights is not None:
        wz = zscores * weights
        csum_wz = np.empty(len(wz) + 1, dtype=np.float64)
        csum_wz[0] = 0
        np.cumsum(wz, out=csum_wz[1:])
        csum_w = np.empty(len(weights) + 1, dtype=np.float64)
        csum_w[0] = 0
        np.cumsum(weights, out=csum_w[1:])
        num = csum_wz[bin_ends] - csum_wz[bin_ends - window_size]
        den = csum_w[bin_ends] - csum_w[bin_ends - window_size]
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(den > 0, num / den, np.nan)
    else:
        csum = np.empty(len(zscores) + 1, dtype=np.float64)
        csum[0] = 0
        np.cumsum(zscores, out=csum[1:])
        return (csum[bin_ends] - csum[bin_ends - window_size]) / window_size
```

**Unweighted case** (point estimates): Build a prefix sum of z-scores. The mean of window `[e - W, e)` is `(csum[e] - csum[e - W]) / W`.

**Weighted case** (bootstrap): Build prefix sums of both `z * w` and `w`. The weighted mean is `sum(z*w) / sum(w)` over the window.

---

## 3. Bin Identities: What Each Bin Represents

Each bin (window position) is identified by its `bin_end` value — the 1-indexed right-edge position in the rank-sorted variant list. This means:

| `bin_end` | Variants included | Interpretation |
|-----------|-------------------|----------------|
| 1000 | Ranks 1–1000 | Top 1000 most deleterious |
| 1010 | Ranks 11–1010 | Shifted by 10 |
| 1020 | Ranks 21–1020 | Shifted by 20 |
| ... | ... | ... |
| N | Ranks N-999 to N | Least deleterious window |

The x-axis of the final plot is `log10(1 / bin_end)`, so:
- The **left** side of the plot (large negative x) corresponds to large `bin_end` values = windows deep into the low-pathogenicity tail.
- The **right** side (small negative x) corresponds to small `bin_end` values = windows among the top-ranked, most deleterious variants.

---

## 4. Bootstrap Resampling: How It Works

The bootstrap resamples **genes (regions)**, not individual variants. This accounts for the fact that variants within the same gene are not independent.

### Setup: Region Index Mapping

Each gene is mapped to an integer index for efficient array-based resampling:

```python
all_regions = sorted(gene_trait_df['region'].unique().to_list())
n_regions = len(all_regions)
region_idx_map = pl.DataFrame({
    'region': all_regions,
    '_region_idx': np.arange(n_regions, dtype=np.int32)
})
ranked_zscores = ranked_zscores.join(region_idx_map, on='region', how='left')
```

Each variant in `ranked_zscores` now carries its gene's integer index (`_region_idx`).

### Bootstrap Loop

For each of 1000 bootstrap iterations:

1. **Resample genes with replacement**: Draw `n_regions` gene indices from `[0, n_regions)` with replacement.
2. **Compute region weights**: Use `np.bincount` to count how many times each gene was selected. This produces a weight vector of length `n_regions` — a gene drawn twice gets weight 2, a gene not drawn gets weight 0.
3. **Map to variant weights**: Each variant inherits the weight of its gene via `region_weights[d['region_idx']]`.
4. **Compute weighted sliding window means**: Call `sliding_window_means` with the per-variant weights, which computes `sum(z * w) / sum(w)` over each window.

```python
rng = np.random.default_rng(42)

for b in tqdm(range(n_boot), desc='Bootstrapping'):
    idx = rng.choice(n_regions, size=n_regions, replace=True)
    region_weights = np.bincount(idx, minlength=n_regions).astype(np.float64)

    for annotation in annotations_list:
        d = anno_arrays[annotation]
        variant_weights = region_weights[d['region_idx']]
        boot_means[annotation][b] = sliding_window_means(
            d['zscores'], d['bin_ends'], window_size, weights=variant_weights
        )
```

This is equivalent to a **cluster bootstrap** where genes are the clusters, but implemented efficiently via weighted means rather than literal subsetting.

---

## 5. Confidence Intervals: How They Are Computed

After all 1000 bootstrap iterations, the notebook computes **percentile-based 95% confidence intervals** and the bootstrap standard error at each window position:

```python
all_results.append(pl.DataFrame({
    'annotation': annotation,
    'bin_end': d['bin_ends'].astype(np.float64),
    'mean_zscore': point_means[annotation].astype(np.float32),
    'se_zscore': np.nanstd(bm, axis=0).astype(np.float32),
    'ci_lower': np.nanpercentile(bm, 2.5, axis=0).astype(np.float32),
    'ci_upper': np.nanpercentile(bm, 97.5, axis=0).astype(np.float32),
}))
```

| Statistic | Computation | Meaning |
|-----------|-------------|---------|
| `mean_zscore` | Unweighted sliding window mean (point estimate) | Best estimate of mean z-score in the window |
| `se_zscore` | `np.nanstd(boot_means, axis=0)` | Standard deviation of the 1000 bootstrap means |
| `ci_lower` | `np.nanpercentile(boot_means, 2.5, axis=0)` | 2.5th percentile of bootstrap distribution |
| `ci_upper` | `np.nanpercentile(boot_means, 97.5, axis=0)` | 97.5th percentile of bootstrap distribution |

The confidence intervals are plotted as semi-transparent ribbons (`geom_ribbon`) around the point estimate lines.

---

## Summary of the Full Pipeline

```
Gene-trait associations (670 genes)
        │
        ├──→ Annotation scores (missense SNPs in LIP domains)
        │       │
        │       ├── Direction-correct scores
        │       ├── Remove imputed (fill-NA) scores (anti-join)
        │       └── Rank variants per annotation (most deleterious = rank 1)
        │
        ├──→ Per-variant phenotype z-scores (MAC ≤ 20)
        │       │
        │       └── Direction-correct by LOFTEE correlation sign
        │
        └──→ Join z-scores ↔ annotation ranks on (id, region)
                │
                ├── Sort by rank per annotation
                ├── Compute unweighted sliding window means (prefix sums)
                ├── Bootstrap: resample genes 1000×, compute weighted window means
                └── Output: mean_zscore, se, ci_lower, ci_upper per (annotation, bin_end)
```
