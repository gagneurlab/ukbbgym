# Odds Ratio Computation with Asof Join

This document explains the algorithm used in `pheno_odds_across_genes.ipynb` to compute odds ratios efficiently using cumulative sums and asof joins instead of cross joins.

## The Goal

For each (annotation, rank_cutoff), we want to build a 2×2 contingency table:

```
                    | extreme pheno | not extreme
-------------------------------------------------
rank <= cutoff      |      a        |      b
rank > cutoff       |      c        |      d
```

Then compute: **OR = (a/b) / (c/d)**

Where:
- **rank**: variant's annotation score rank (1 = highest score)
- **extreme pheno**: variant associated with phenotype in top 1%
- **cutoff**: we test 101 different rank cutoffs (e.g., top 1000, top 2000, ..., top 10M variants)

## The Naive Approach (Cross Join)

```python
# Replicate every variant 101 times (one per cutoff)
gp.join(cutoffs, how='cross')
  .group_by(['annotation', 'rank_cutoff'])
  .agg(
      a = ((rank <= cutoff) & is_extreme).sum(),
      b = ((rank <= cutoff) & ~is_extreme).sum(),
      ...
  )
```

**Problem:** With ~18M variants × 11 annotations × 101 cutoffs = **20 billion rows** in the intermediate result. Even with streaming, this is slow and memory-intensive.

## The Efficient Approach (Cumulative Sums + Asof Join)

**Key insight:** Since we're testing rank cutoffs in sorted order, we can use **cumulative sums** to avoid replicating data.

---

## Step-by-Step Walkthrough

### Step 1: Build Base Data

```python
gp_base = (
    appv
    .join(id_region, on='id', how='inner')       # add gene region to each variant
    .join(gene_trait_df.lazy(), on=['region', 'phenotype'], how='inner')
                                                  # keep only matching gene-phenotype pairs
    .with_columns(is_extreme = (...))             # True if phenotype in top 1%
    .select(['id', 'region', 'is_extreme'])
    .collect(engine='streaming')
)
```

**Result:** ~18M rows. Each variant has a boolean `is_extreme` flag.

**Example:**
```
id              | region          | is_extreme
chr2:21012073   | ENSG00000084674 | False
chr17:1474584   | ENSG00000197879 | True
...
```

---

### Step 2: Extract One Annotation's Ranks

We process annotations one at a time to control memory usage.

```python
anno_ranks = (
    melted_anno
    .filter(pl.col('annotation') == 'am_pathogenicity')
    .select(['id', 'region', 'annotation_score_dircor_rank_desc'])
)
```

**Example:**
```
id              | region          | rank
chr2:21012073   | ENSG00000084674 | 535
chr17:1474584   | ENSG00000197879 | 17576
...
```

---

### Step 3: Aggregate to Unique Ranks

```python
rank_counts = (
    gp_base.lazy()
    .join(anno_ranks.lazy(), on=['id', 'region'], how='inner')
    .group_by('annotation_score_dircor_rank_desc')
    .agg(
        n_extreme = pl.col('is_extreme').sum().cast(pl.Int64),
        n_total = pl.len().cast(pl.Int64),
    )
    .collect()
    .sort('annotation_score_dircor_rank_desc')
)
```

Variants with the same rank (ties from `method="max"`) are collapsed into counts.

**Example:**
```
rank | n_extreme | n_total
  1  |     0     |    1       ← 1 variant at rank 1, not extreme
  2  |     1     |    1       ← 1 variant at rank 2, is extreme
  3  |     0     |    3       ← 3 tied variants at rank 3, none extreme
  7  |     1     |    2       ← 2 variants tied at rank 7, one extreme
  ...
```

---

### Step 4: Compute Cumulative Sums (The Key Trick!)

```python
.with_columns(
    cum_extreme = pl.col('n_extreme').cum_sum(),
    cum_total = pl.col('n_total').cum_sum(),
)
```

After sorting by rank ascending, cumulative sums give us running totals:

```
rank | n_extreme | n_total | cum_extreme | cum_total
  1  |     0     |    1    |      0      |     1       ← "top 1" has 0 extreme, 1 total
  2  |     1     |    1    |      1      |     2       ← "top 2" has 1 extreme, 2 total
  3  |     0     |    3    |      1      |     5       ← "top 5" has 1 extreme, 5 total
  7  |     1     |    2    |      2      |     7       ← "top 7" has 2 extreme, 7 total
  ...
```

**Critical insight:**
- `cum_extreme` at rank R = number of extreme variants with rank ≤ R = **cell `a`** in the 2×2 table
- `cum_total` at rank R = total variants with rank ≤ R = **`a + b`** in the 2×2 table

---

### Step 5: Store Grand Totals

```python
total_extreme = rank_counts['n_extreme'].sum()   # a + c (all extreme variants)
total_count = rank_counts['n_total'].sum()        # a + b + c + d (all variants)
```

These are needed to compute cells `c` and `d` by subtraction.

---

### Step 6: Asof Join — Look Up Cumulative Values at Each Cutoff

```python
result = (
    cutoffs_sorted
    .join_asof(
        rank_counts,
        left_on='rank_cutoff',
        right_on='annotation_score_dircor_rank_desc',
        strategy='backward',
    )
)
```

#### How `join_asof` Works

For each row in the **left** table (cutoffs), it finds the nearest matching row in the **right** table (rank_counts) where `right_key <= left_key`. Both tables must be sorted on their join keys.

**Example** with cutoffs [1000, 2000, 5000] and ranks [..., 998, 999, 1003, 1005, ...]:

```
cutoff=1000 → finds rank=999 (last rank ≤ 1000) → gets cum_extreme=X, cum_total=Y at rank 999
cutoff=2000 → finds rank=1005 (last rank ≤ 2000) → gets cumulative counts at rank 1005
cutoff=5000 → finds rank=4999 (last rank ≤ 5000) → gets cumulative counts at rank 4999
```

This is essentially a **binary search** for each cutoff — O(log N) per lookup instead of scanning all rows.

**Why `strategy='backward'`:** We want the **last** rank that is ≤ cutoff (looking backward), not the next one after it.

**Why `fill_null(0)`:** If cutoff is smaller than the minimum rank, no match is found → null. This means zero variants are "above the cutoff," so we fill with 0.

```python
.with_columns(
    cum_extreme = pl.col('cum_extreme').fill_null(0),
    cum_total = pl.col('cum_total').fill_null(0),
)
```

---

### Step 7: Compute 2×2 Table Cells

```python
.with_columns(
    n_dis_above    = pl.col('cum_extreme'),                           # a
    n_notdis_above = pl.col('cum_total') - pl.col('cum_extreme'),     # b = (a+b) - a
    n_dis_below    = pl.lit(total_extreme) - pl.col('cum_extreme'),   # c = (a+c) - a
    n_notdis_below = pl.lit(total_count - total_extreme)              # d = (b+d) - b
                   - (pl.col('cum_total') - pl.col('cum_extreme')),
)
```

All four cells are derived from:
- `cum_extreme` (from asof join)
- `cum_total` (from asof join)
- `total_extreme` (grand total)
- `total_count` (grand total)

**Example** for cutoff=1000 where `cum_extreme=10`, `cum_total=950`, `total_extreme=200`, `total_count=18000000`:

```
a = 10                          (10 extreme variants in top 1000)
b = 950 - 10 = 940              (940 non-extreme variants in top 1000)
c = 200 - 10 = 190              (190 extreme variants below top 1000)
d = (18000000-200) - (950-10)   (rest are non-extreme below top 1000)
  = 17999800 - 940 = 17998860
```

---

### Step 8: Combine Results and Compute Odds Ratios

After looping through all 11 annotations, we have 11 DataFrames of 101 rows each (1,111 total rows).

```python
or_df = (
    pl.concat(results)
    .with_columns(
        odds_ratio = (pl.col('n_dis_above') / pl.col('n_notdis_above'))
                   / (pl.col('n_dis_below') / pl.col('n_notdis_below'))
    )
    .filter(pl.col('odds_ratio').is_finite() & (pl.col('odds_ratio') > 0))
    .with_columns(
        se_log_or = (
            1/pl.col('n_dis_above').cast(pl.Float64)
            + 1/pl.col('n_notdis_above').cast(pl.Float64)
            + 1/pl.col('n_dis_below').cast(pl.Float64)
            + 1/pl.col('n_notdis_below').cast(pl.Float64)
        ).sqrt(),
    )
    .with_columns(
        ci_lower = (pl.col('odds_ratio').log() - 1.96 * pl.col('se_log_or')).exp(),
        ci_upper = (pl.col('odds_ratio').log() + 1.96 * pl.col('se_log_or')).exp(),
    )
)
```

**Confidence Intervals:** We use the **Woolf method** for log-OR standard error:

```
SE(ln OR) = sqrt(1/a + 1/b + 1/c + 1/d)
95% CI = exp(ln(OR) ± 1.96 * SE)
```

We drop rows where OR is infinite or zero (degenerate cases with zero cells).

---

## Performance Comparison

| Approach | Intermediate Rows | Memory | Speed |
|----------|------------------|--------|-------|
| **Cross join** | 18M × 11 × 101 = **20 billion** | Very high | Slow even with streaming |
| **Asof join** | 18M per annotation (processed sequentially) | Moderate (~1-1.5GB peak) | **~10x faster** |

---

## Why Process One Annotation at a Time?

With 11 annotations × ~18M variants, storing all `rank_counts` at once would require ~200M rows in memory (~8GB). By looping, we only hold one annotation's data (~18M rows, ~1.5GB) at a time, with each intermediate DataFrame freed after computing that annotation's results.

The final result is only 1,111 rows (11 annotations × 101 cutoffs), which is tiny.

---

## Summary

Instead of replicating each of ~18M rows 101 times (cross join → 1.8B rows per annotation), we:
1. Sort variants by rank once
2. Compute cumulative sums once
3. Do 101 binary search lookups (asof join)

**Same result, ~1000x less work.**
