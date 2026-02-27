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

#### The problem

A naive approach would loop over every window and sum the z-scores inside it:

```
for each window position:
    mean = sum(zscores[start:end]) / window_size   # O(window_size) per window
```

With ~30,000+ window positions and `window_size=1000`, that is ~30 million additions. The prefix sum trick reduces each window to a single subtraction.

#### What is a prefix sum?

A **prefix sum** (or cumulative sum) array `csum` is defined so that `csum[i]` = sum of the first `i` elements:

```
zscores:  [ z0,   z1,   z2,   z3,   z4,  ... ]
csum:  [ 0,  z0,  z0+z1, z0+z1+z2, z0+z1+z2+z3, ... ]
         ^
         csum[0] = 0 (sentinel)
```

The key identity: **the sum of any contiguous slice `[a, b)` is just `csum[b] - csum[a]`**. No matter how wide the window, computing its sum is always one subtraction — O(1).

#### Concrete walkthrough (unweighted)

Suppose `window_size = 3` and we have 6 z-scores sorted by annotation rank:

```
index:     0     1     2     3     4     5
zscores: [0.5,  0.3,  0.4,  0.1,  0.2, -0.1]
csum:    [0.0,  0.5,  0.8,  1.2,  1.3,  1.5,  1.4]
               ^csum[1]                      ^csum[6]
```

`bin_ends = [3, 4, 5, 6]` (with `step_size=1` here for illustration).

| bin_end | window slice | sum via prefix sums | mean |
|---------|-------------|---------------------|------|
| 3 | `zscores[0:3]` = [0.5, 0.3, 0.4] | `csum[3] - csum[0]` = 1.2 - 0.0 = 1.2 | 0.40 |
| 4 | `zscores[1:4]` = [0.3, 0.4, 0.1] | `csum[4] - csum[1]` = 1.3 - 0.5 = 0.8 | 0.27 |
| 5 | `zscores[2:5]` = [0.4, 0.1, 0.2] | `csum[5] - csum[2]` = 1.5 - 0.8 = 0.7 | 0.23 |
| 6 | `zscores[3:6]` = [0.1, 0.2, -0.1] | `csum[6] - csum[3]` = 1.4 - 1.2 = 0.2 | 0.07 |

Each window mean costs one subtraction and one division, regardless of `window_size`.

#### The code

```python
def sliding_window_means(zscores, bin_ends, window_size, weights=None):
    if weights is not None:
        # --- Weighted case (used during bootstrap) ---
        wz = zscores * weights                          # element-wise z * w
        csum_wz = np.empty(len(wz) + 1, dtype=np.float64)
        csum_wz[0] = 0
        np.cumsum(wz, out=csum_wz[1:])                  # prefix sum of (z * w)
        csum_w = np.empty(len(weights) + 1, dtype=np.float64)
        csum_w[0] = 0
        np.cumsum(weights, out=csum_w[1:])               # prefix sum of w
        num = csum_wz[bin_ends] - csum_wz[bin_ends - window_size]  # sum(z*w) in window
        den = csum_w[bin_ends] - csum_w[bin_ends - window_size]    # sum(w) in window
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(den > 0, num / den, np.nan)
    else:
        # --- Unweighted case (used for point estimates) ---
        csum = np.empty(len(zscores) + 1, dtype=np.float64)
        csum[0] = 0
        np.cumsum(zscores, out=csum[1:])                 # prefix sum of z
        return (csum[bin_ends] - csum[bin_ends - window_size]) / window_size
```

**Unweighted case** (point estimates): Build one prefix sum of z-scores. The mean of window `[e - W, e)` is `(csum[e] - csum[e - W]) / W`.

**Weighted case** (bootstrap): Build two prefix sums — one for `z * w` and one for `w`. The weighted mean of window `[e - W, e)` is `sum(z*w) / sum(w)` = `(csum_wz[e] - csum_wz[e - W]) / (csum_w[e] - csum_w[e - W])`. If the denominator is zero (all variants in the window have weight 0), the result is `NaN`.

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

The bootstrap resamples **genes (regions)**, not individual variants. This accounts for the fact that variants within the same gene are not independent (they share the same gene-trait association, the same set of carriers, etc.).

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

Each variant in `ranked_zscores` now carries its gene's integer index (`_region_idx`). This allows the bootstrap loop to work entirely with fast numpy integer arrays rather than string comparisons.

### Bootstrap Loop

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

For each of 1000 bootstrap iterations:

1. **Resample genes with replacement**: Draw `n_regions` gene indices from `[0, n_regions)` with replacement.
2. **Compute region weights**: `np.bincount(idx, minlength=n_regions)` counts how many times each gene was drawn. A gene drawn twice gets weight 2, a gene not drawn gets weight 0.
3. **Map to variant weights**: Each variant inherits the weight of its parent gene: `variant_weights = region_weights[d['region_idx']]`. So if gene G was sampled 3 times, every variant from gene G gets weight 3.
4. **Compute weighted sliding window means**: The `sliding_window_means` function computes `sum(z * w) / sum(w)` over each window using the prefix sum trick (see section 2).

### Reweighting vs. Re-ranking: Why the Rank Order Is Fixed

You might expect that a bootstrap resample should produce a **new ranking** — after all, if some genes are dropped (weight 0) and others duplicated (weight 2+), shouldn't the variant positions shift? Here is why the notebook does **not** re-rank, and what that design choice implies.

#### What a literal "re-rank" bootstrap would do

1. Resample genes with replacement.
2. Take all variants from the sampled genes (with duplicates for genes sampled >1 time).
3. **Re-sort** these variants by annotation score and assign new ranks.
4. Compute sliding windows over the new ranking.

Because the annotation scores are intrinsic to each variant (not gene-dependent), duplicating a gene's variants just inserts tied copies at the same score. Removing a gene's variants creates gaps. After re-sorting, the relative order of the remaining variants is unchanged — only the positions shift. This is computationally expensive (re-sorting ~300k variants × 1000 iterations × 6 annotations), but conceptually clean.

#### What the notebook actually does (fixed-rank reweighting)

Instead of re-ranking, the notebook:

1. Keeps the original sorted array of variants fixed in its original rank order.
2. Assigns each variant a weight (0, 1, 2, ...) based on how many times its gene was sampled.
3. Computes a **weighted mean** over the same fixed windows.

A variant from a gene that was not sampled (weight = 0) contributes nothing to either the numerator `sum(z*w)` or denominator `sum(w)` — it is effectively invisible. A variant from a gene sampled twice contributes double.

#### The key difference and what it means

The difference shows up at window boundaries. Consider a concrete example:

```
Original rank order (window_size=3):
Position:  1     2     3  |  4     5     6
Gene:      A     A     B  |  B     C     C
z-score:   0.5   0.4   0.3  0.2   0.1   0.0

Window [1-3]: variants from genes A, A, B
Window [4-6]: variants from genes B, C, C
```

Now suppose bootstrap samples genes {A, C} (gene B has weight 0):

**Fixed-rank reweighting** (what the code does):
```
Position:  1     2     3  |  4     5     6
Weight:    1     1     0  |  0     1     1
Window [1-3]: weighted mean = (0.5×1 + 0.4×1 + 0.3×0) / (1+1+0) = 0.45
Window [4-6]: weighted mean = (0.2×0 + 0.1×1 + 0.0×1) / (0+1+1) = 0.05
```

**Re-ranking** (the alternative approach):
```
Remaining variants after removing gene B:
Position:  1     2     3     4
Gene:      A     A     C     C
z-score:   0.5   0.4   0.1   0.0

Window [1-3]: mean of [0.5, 0.4, 0.1] = 0.33
(Gene C's variants have "slid up" into the top window)
```

The fixed-rank approach gives 0.45 for window [1-3]; the re-ranking approach gives 0.33 because gene C's low-z variants have filled the gap left by gene B. In other words:

- **Fixed-rank reweighting** answers: *"At this rank position in the annotation, what is the mean z-score accounting for gene-level sampling uncertainty?"*
- **Re-ranking** answers: *"If the population of genes were different, what would the top-N variants' mean z-score be?"*

#### Why fixed-rank reweighting is the appropriate choice here

The notebook's goal is to estimate **how mean phenotypic effect varies as a function of annotation rank**, with uncertainty bands. The x-axis is "rank position in the annotation," and the question at each x-position is: "how uncertain is the mean z-score at this rank, given that genes are the independent units?"

Fixed-rank reweighting is the standard approach for this because:

1. **The x-axis should remain stable across bootstrap iterations.** Each window position corresponds to a fixed annotation score range. Re-ranking would shift what scores each window corresponds to in every iteration, making the x-axis meaning inconsistent.

2. **It is a standard cluster bootstrap.** This is the textbook "weighted bootstrap" or "Bayesian bootstrap" approach for clustered data (see Davison & Hinkley, 1997). Rather than literally duplicating and removing observations, you assign multinomial weights to clusters. It is mathematically equivalent to the literal approach when computing means and other smooth statistics at fixed positions.

3. **Annotation rank is a variant property, not a gene property.** The annotation score (and therefore the rank) is determined by the variant's protein-level impact. It does not change when you resample genes. Re-ranking would conflate two sources of variation: gene sampling and rank assignment.

#### Potential edge-case to be aware of

If a window happens to contain variants from very few genes, and the bootstrap assigns weight 0 to most of those genes, the denominator `sum(w)` can become very small, producing noisy or `NaN` estimates for that window in that iteration. The code handles the `sum(w) = 0` case explicitly:

```python
with np.errstate(divide='ignore', invalid='ignore'):
    return np.where(den > 0, num / den, np.nan)
```

These `NaN` values are handled downstream by `np.nanstd` and `np.nanpercentile`, which skip them. In practice, with `window_size=1000` spanning variants from many genes, this is rarely an issue except possibly at the extreme tails.

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
