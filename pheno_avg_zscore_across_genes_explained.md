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

#### The code (point estimates)

For the unweighted point estimates, the standard prefix sum approach is used:

```python
def sliding_window_means_unweighted(zscores, bin_ends, window_size):
    csum = np.empty(len(zscores) + 1, dtype=np.float64)
    csum[0] = 0.0
    np.cumsum(zscores, out=csum[1:])
    return (csum[bin_ends] - csum[bin_ends - window_size]) / window_size
```

Build one prefix sum of z-scores. The mean of window `[e - W, e)` is `(csum[e] - csum[e - W]) / W`.

The bootstrap uses a different approach — a **virtual expanded array** with binary search — described in section 4.

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

## 4. Bootstrap Resampling: Re-ranking Bootstrap

The bootstrap resamples **genes (regions)**, not individual variants. This accounts for the fact that variants within the same gene are not independent (they share the same gene-trait association, the same set of carriers, etc.).

Crucially, the notebook uses a **re-ranking bootstrap** rather than a fixed-rank reweighting approach. This means that when genes are resampled, variants from dropped genes are removed and variants from duplicated genes fill their place, causing the rank positions to shift. This correctly captures the uncertainty that arises from gene-dependent annotation score distributions.

### Why re-ranking rather than fixed-rank reweighting?

Annotation scores are **gene-dependent** — a highly conserved gene will have systematically higher pathogenicity scores than a loosely constrained gene. When you resample genes, you change the composition of variants at each rank position. A fixed-rank approach would keep dropped genes' variants "reserving" their rank positions (contributing zero weight), preventing lower-ranked variants from sliding up. This underestimates uncertainty, particularly at the top of the ranking where a few large constrained genes can dominate.

Consider a concrete example:

```
Original rank order (window_size=3):
Position:  1     2     3  |  4     5     6
Gene:      A     A     B  |  B     C     C
z-score:   0.5   0.4   0.3  0.2   0.1   0.0
```

Bootstrap samples genes {A, C} (gene B dropped):

**Fixed-rank reweighting** (keeps rank positions fixed):
```
Position:  1     2     3  |  4     5     6
Weight:    1     1     0  |  0     1     1
Window [1-3]: weighted mean = (0.5×1 + 0.4×1 + 0.3×0) / (1+1+0) = 0.45
```

**Re-ranking** (what the code does — variants slide up to fill gaps):
```
Expanded array after removing gene B:
Position:  1     2     3     4
Gene:      A     A     C     C
z-score:   0.5   0.4   0.1   0.0

Window [1-3]: mean of [0.5, 0.4, 0.1] = 0.33
```

The re-ranking approach gives 0.33 because gene C's lower-z variants fill the gap left by gene B. This correctly reflects: *"if gene B weren't in our dataset, the top-3 window would include gene C's variants."* The fixed-rank approach gives 0.45, which overstates the signal by ignoring that removing gene B changes which variants constitute the "top 3."

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

### The Virtual Expanded Array (converting time to memory)

The naive re-ranking approach would be: for each bootstrap iteration, build a new array by repeating/removing variants according to gene weights, then re-sort and apply the sliding window. Re-sorting ~300k variants × 1000 iterations × 6 annotations would be prohibitively slow.

The key insight is that **re-sorting is unnecessary**. Because annotation scores are intrinsic to each variant, the relative order never changes — only the positions shift. When gene B is dropped, gene C's variants don't change their scores; they just move up to fill the gap. When gene A is duplicated, its variants appear twice at the same positions in the sort order.

This means the "re-ranked" array is equivalent to the original sorted array with each variant repeated `w[i]` times (where `w[i]` is the bootstrap weight of variant `i`'s gene). We call this the **virtual expanded array** — we never materialize it, but compute its prefix sum analytically using binary search.

#### How it works

For each bootstrap iteration, we compute two cumulative arrays from the original sorted variants:

```python
# vw[i] = bootstrap weight of variant i's gene (0, 1, 2, ...)
vw = region_weights[d['region_idx']]

# cw[i] = cumulative weight up to variant i  (cw[0] = 0)
# This tells us: variant i's copies span expanded positions (cw[i-1]+1) to cw[i]
cw[0] = 0.0
np.cumsum(vw, out=cw[1:])

# cwz[i] = cumulative weighted z-score up to variant i  (cwz[0] = 0)
cwz[0] = 0.0
np.cumsum(z * vw, out=cwz[1:])
```

The total length of the virtual expanded array is `cw[N]`.

To compute the prefix sum of the expanded array at any position `e`, we use binary search:

```python
def expanded_csum_at(positions, cw, cwz, z, N):
    # Find which original variant spans each expanded position
    k = np.searchsorted(cw, positions, side='left')
    k_safe = np.clip(k, 1, N)
    # Full prefix sum up to variant k-2, plus partial contribution from variant k-1
    result = cwz[k_safe - 1] + z[k_safe - 1] * (positions - cw[k_safe - 1])
    result = np.where(positions <= 0, 0.0, result)
    return result
```

#### Concrete walkthrough

```
Original sorted array:
Index:    0     1     2     3
z-score:  0.5   0.4   0.3   0.1
Gene:     A     A     B     B
Weight:   1     1     0     0    (bootstrap sampled {A, C} but no C variants here)

cw  = [0, 1, 2, 2, 2]     ← variant 2 and 3 have weight 0, so cw stays flat
cwz = [0, 0.5, 0.9, 0.9, 0.9]

Virtual expanded array: [0.5, 0.4]  (only 2 elements, total_expanded = 2)
```

Query `expanded_csum_at(e=1)`:
- `searchsorted(cw, 1, 'left')` → k=1 (cw[1]=1 ≥ 1)
- result = cwz[0] + z[0] × (1 - cw[0]) = 0 + 0.5 × 1 = 0.5 ✓

Query `expanded_csum_at(e=2)`:
- `searchsorted(cw, 2, 'left')` → k=2 (cw[2]=2 ≥ 2)
- result = cwz[1] + z[1] × (2 - cw[1]) = 0.5 + 0.4 × 1 = 0.9 ✓

With partial copies (gene A sampled twice, weight=2):
```
Weight:   2     2     0     0
cw  = [0, 2, 4, 4, 4]
cwz = [0, 1.0, 1.8, 1.8, 1.8]

Virtual expanded array: [0.5, 0.5, 0.4, 0.4]  (total_expanded = 4)
```

Query `expanded_csum_at(e=3)`:
- `searchsorted([0,2,4,4,4], 3, 'left')` → k=2 (cw[2]=4 ≥ 3)
- result = cwz[1] + z[1] × (3 - cw[1]) = 1.0 + 0.4 × (3-2) = 1.4
- Expected: 0.5 + 0.5 + 0.4 = 1.4 ✓

The `searchsorted` handles the "partial copy" case — when the window boundary falls in the middle of a variant's repeated copies.

#### The sliding window on the expanded array

```python
def sliding_window_reranked(z, cw, cwz, N, total_expanded, bin_ends, window_size):
    valid_mask = bin_ends <= total_expanded
    valid_bins = bin_ends[valid_mask]
    means = np.full(len(bin_ends), np.nan, dtype=np.float64)
    if len(valid_bins) == 0:
        return means
    right = expanded_csum_at(valid_bins, cw, cwz, z, N)
    left = expanded_csum_at(valid_bins - window_size, cw, cwz, z, N)
    means[valid_mask] = (right - left) / window_size
    return means
```

The window mean at expanded position `e` is `(expanded_csum[e] - expanded_csum[e - W]) / W`. Bin positions that exceed the expanded array length for a given bootstrap iteration are set to `NaN`.

### Performance: Why This Is Fast

The approach avoids both re-sorting and materializing the expanded array:

| Operation | Cost |
|-----------|------|
| `region_weights[region_idx]` — map gene weights to variants | O(N) |
| `cumsum(vw)` and `cumsum(z*vw)` — build cw and cwz | O(N) |
| `searchsorted` — binary search for ~30k bin positions | O(n_bins × log N) |
| **Total per annotation per iteration** | **O(N + n_bins × log N)** |

Pre-allocated buffers (`cw_buf`, `cwz_buf`) avoid memory allocation in the hot loop:

```python
max_N = max(d['N'] for d in anno_arrays.values())
cw_buf = np.empty(max_N + 1, dtype=np.float64)
cwz_buf = np.empty(max_N + 1, dtype=np.float64)
```

### Bootstrap Loop

```python
rng = np.random.default_rng(42)

for b in tqdm(range(n_boot), desc='Re-ranking bootstrap'):
    idx = rng.choice(n_regions, size=n_regions, replace=True)
    region_weights = np.bincount(idx, minlength=n_regions).astype(np.float64)

    for annotation in annotations_list:
        d = anno_arrays[annotation]
        N = d['N']
        z = d['zscores']
        vw = region_weights[d['region_idx']]

        cw = cw_buf[:N + 1]
        cw[0] = 0.0
        np.cumsum(vw, out=cw[1:])

        cwz = cwz_buf[:N + 1]
        cwz[0] = 0.0
        np.cumsum(z * vw, out=cwz[1:])

        total_expanded = int(cw[N])

        boot_means[annotation][b] = sliding_window_reranked(
            z, cw, cwz, N, total_expanded, d['bin_ends'], window_size
        )
```

For each of 1000 bootstrap iterations:

1. **Resample genes with replacement**: Draw `n_regions` gene indices with replacement.
2. **Compute region weights**: `np.bincount` counts how many times each gene was drawn.
3. **Map to variant weights**: Each variant inherits the weight of its parent gene.
4. **Build cumulative arrays**: `cw` (cumulative weights) and `cwz` (cumulative weighted z-scores) using pre-allocated buffers.
5. **Compute sliding window means**: Binary search into `cw` to evaluate the expanded prefix sum at each window edge, without materializing the expanded array.

---

## 5. Confidence Intervals: How They Are Computed

After all 1000 bootstrap iterations, the notebook computes **percentile-based 95% confidence intervals**:

```python
pm = point_means[annotation].astype(np.float32)

all_results.append(pl.DataFrame({
    'annotation': annotation,
    'bin_end': d['bin_ends'].astype(np.float64),
    'mean_zscore': pm,
    'ci_lower': np.nanpercentile(bm, 2.5, axis=0).astype(np.float32),
    'ci_upper': np.nanpercentile(bm, 97.5, axis=0).astype(np.float32),
}))
```

| Statistic | Computation | Meaning |
|-----------|-------------|---------|
| `mean_zscore` | Unweighted sliding window mean (point estimate on full data) | Best estimate of mean z-score in the window |
| `ci_lower` | `np.nanpercentile(boot_means, 2.5, axis=0)` | 2.5th percentile of 1000 bootstrap means |
| `ci_upper` | `np.nanpercentile(boot_means, 97.5, axis=0)` | 97.5th percentile of 1000 bootstrap means |

### Why percentile-based CIs?

Percentile-based CIs are a natural choice for the bootstrap: they directly reflect the distribution of the bootstrap replicates, making no assumptions about normality or symmetry. If the bootstrap distribution is skewed at certain rank positions (e.g. at the top of the ranking where a few large genes dominate), the percentile CI will correctly capture that asymmetry.

### NaN handling at the tail

With the re-ranking bootstrap, the virtual expanded array has a different length in each iteration (since some genes are dropped and others duplicated). When the expanded array is shorter than the original variant count N, bin positions near the tail exceed the expanded array length and are set to `NaN` for that iteration. However, every bin that fits within the expanded array contains exactly `window_size` variants — there are no sparse or partially-filled bins. `np.nanpercentile` simply ignores the `NaN` values, computing the percentile over whichever iterations had enough expanded variants to reach that bin position.

The confidence intervals are plotted as semi-transparent ribbons (`geom_ribbon`) around the point estimate lines.

---

## Summary of the Full Pipeline

```
Gene-trait associations (670 genes)
        │
        ├──→ Annotation scores (missense SNPs)
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
                ├── Re-ranking bootstrap: resample genes 1000×
                │       │
                │       ├── Build virtual expanded array (cumulative weights + binary search)
                │       └── Sliding window on expanded array (no materialization)
                │
                └── Output: mean_zscore, ci_lower, ci_upper per (annotation, bin_end)
                            (Percentile CIs: 2.5th / 97.5th of bootstrap means)
```
