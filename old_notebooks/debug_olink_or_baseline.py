#!/usr/bin/env python3
"""
Debug script to understand why olink continuous OR doesn't reach category OR
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the annotation data (missense-filtered)
print("Loading annotation data...")
anno = pl.read_parquet('PATH_TO_FILE')

# Load olink whitelist
burden_test = pl.read_parquet('PATH_TO_FILE')
olink_whitelist = burden_test.rename({'gene': 'region'}).select(['region']).unique()

print(f"Olink whitelist: {olink_whitelist.shape[0]} genes")

# Filter to missense SNPs in olink genes (same as notebook)
anno_missense = (
    anno
    .filter(
        (pl.col('region').is_in(olink_whitelist['region'])) &
        (pl.col('consequence_missense_variant') == True) &
        (pl.col('ref').str.len_chars() == 1) &
        (pl.col('alt').str.len_chars() == 1)
    )
    .select(['id', 'region', 'cadd_raw'])
    .drop_nulls('cadd_raw')
)

print(f"\nMissense SNPs: {anno_missense.shape[0]:,}")

# Load APPV
print("\nLoading APPV data...")
appv = pl.read_parquet('PATH_TO_FILE')

# Filter APPV (same as notebook)
appv_filtered = appv.filter(pl.col('n_individuals') <= 20)
print(f"APPV filtered (≤20 individuals): {appv_filtered.shape[0]:,}")

# Join with missense annotations
appv_missense = (
    appv_filtered
    .join(anno_missense.select(['id', 'region']), on=['id', 'region'], how='inner')
)

print(f"APPV-missense joined: {appv_missense.shape[0]:,}")

# Compute extreme flag (following notebook logic)
olink_corrs = pl.read_parquet('PATH_TO_FILE')
olink_corrs = olink_corrs.with_columns(
    loftee_corr_dir = pl.col('correlation') / pl.col('correlation').abs()
).select(['region', 'phenotype', 'loftee_corr_dir']).drop_nans()

# Join correlations to determine direction
appv_with_dir = (
    appv_missense
    .join(
        olink_corrs.with_columns(phenotype = pl.col('region') + '_olink'),
        on=['region', 'phenotype'],
        how='left'
    )
    .with_columns(
        # For olink, extreme_pheno_dir = 'bottom', so no flip
        is_extreme = pl.col('mean_pheno_value_ptile') >= 0.99
    )
)

print(f"\n=== BASELINE EXTREME RATES ===")
total = appv_with_dir.shape[0]
extreme = appv_with_dir.filter(pl.col('is_extreme')).shape[0]
print(f"Missense variants extreme rate: {extreme / total * 100:.2f}% ({extreme:,} / {total:,})")

# Now let's see the distribution by CADD score
# Add CADD scores to the data
appv_with_cadd = (
    appv_with_dir
    .join(anno_missense.select(['id', 'region', 'cadd_raw']), on=['id', 'region'], how='left')
    .drop_nulls('cadd_raw')
)

print(f"\nWith CADD scores: {appv_with_cadd.shape[0]:,}")

# Bin by CADD score deciles
appv_with_cadd = appv_with_cadd.with_columns(
    cadd_decile = pl.col('cadd_raw').qcut(10, labels=[f"D{i}" for i in range(1, 11)])
)

# Compute extreme rate per decile
extreme_by_decile = (
    appv_with_cadd
    .group_by('cadd_decile')
    .agg(
        n_total = pl.len(),
        n_extreme = pl.col('is_extreme').sum(),
        cadd_min = pl.col('cadd_raw').min(),
        cadd_max = pl.col('cadd_raw').max(),
    )
    .with_columns(
        extreme_rate = (pl.col('n_extreme') / pl.col('n_total') * 100)
    )
    .sort('cadd_decile')
)

print("\n=== EXTREME RATE BY CADD DECILE ===")
print(extreme_by_decile)

# Create visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Extreme rate by CADD decile
ax1 = axes[0]
decile_data = extreme_by_decile.to_pandas()
ax1.bar(range(len(decile_data)), decile_data['extreme_rate'], color='steelblue')
ax1.axhline(y=1.0, color='red', linestyle='--', label='Baseline (1%)')
ax1.axhline(y=extreme / total * 100, color='orange', linestyle='--',
            label=f'Overall missense ({extreme / total * 100:.1f}%)')
ax1.set_xlabel('CADD Score Decile (1=lowest, 10=highest)', fontsize=12)
ax1.set_ylabel('Extreme Phenotype Rate (%)', fontsize=12)
ax1.set_title('Olink Extreme Rate by CADD Decile\n(Missense variants only)', fontsize=13)
ax1.set_xticks(range(len(decile_data)))
ax1.set_xticklabels([f"D{i}" for i in range(1, 11)])
ax1.legend()
ax1.grid(axis='y', alpha=0.3)

# Plot 2: CADD distribution
ax2 = axes[1]
cadd_values = appv_with_cadd['cadd_raw'].to_numpy()
ax2.hist(cadd_values, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
ax2.set_xlabel('CADD Raw Score', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title('Distribution of CADD Scores\n(Missense variants in olink genes)', fontsize=13)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('PATH_TO_FILE', dpi=150, bbox_inches='tight')
print("\n✓ Saved plot to: PATH_TO_FILE")

# Now compute what the OR should be for continuous approach
print("\n=== COMPUTING EXPECTED OR ===")

# Top 10% CADD vs bottom 10% CADD
top_decile = extreme_by_decile.filter(pl.col('cadd_decile') == 'D10')
bottom_decile = extreme_by_decile.filter(pl.col('cadd_decile') == 'D1')

top_rate = top_decile['extreme_rate'].item()
bottom_rate = bottom_decile['extreme_rate'].item()

print(f"Top 10% CADD extreme rate: {top_rate:.2f}%")
print(f"Bottom 10% CADD extreme rate: {bottom_rate:.2f}%")
print(f"Continuous OR (top vs bottom): {top_rate / bottom_rate:.2f}")
print(f"\nFor comparison:")
print(f"Category OR (missense vs baseline 1%): {(extreme / total * 100) / 1.0:.2f}")
