import os
import gc
import zarr
import yaml
import json
import numba
import numpy as np
import polars as pl
import pandas as pd
from tqdm import tqdm
import statsmodels.api as sm


@numba.njit(parallel=True, fastmath=True)
def _fast_clip_and_sum_allels(arr):
    # Get the shape of the input array
    # Using specific dimensions for clarity with this problem
    n_samples, n_variants, _ = arr.shape
    
    # The sum of two positive int8s can be up to 254. 
    # An int16 is a safe and fast output type.
    output = np.empty((n_samples, n_variants), dtype=np.int8)
    
    # Numba's prange enables automatic parallelization across all your CPU cores
    for i in numba.prange(n_samples):
        for j in range(n_variants):
            # Read two values, perform logic, write one value.
            # This is the "fused" operation.
            val1 = arr[i, j, 0]
            val2 = arr[i, j, 1]
            
            s = 0
            # Since input is int8, this check is faster than max(0, val)
            if val1 > 0:
                s += val1
            if val2 > 0:
                s += val2
            
            output[i, j] = s
            
    return output

def process_genotype_chunk(
    geno: np.array, 
    var_ids: np.array, 
    sample_list: np.array,
    melted_pheno_df: pl.LazyFrame,
    homozygous: bool = False,
    debug: bool = False,
) -> pl.LazyFrame:
    """
    Extract genotypes for a specific gene and return as lazy DataFrame
    """
    geno_clipped = _fast_clip_and_sum_allels(geno)
    # Find heterozygous genotypes (genotype == 1)
    rows, cols = np.where(geno_clipped == 1)
    geno_melt = pl.DataFrame({
        'id': var_ids[rows],
        'individual': sample_list[cols],
        'genotype': 1
    })
    
    # Find homozygous genotypes (genotype == 2)
    if homozygous:
        rows, cols = np.where(geno_clipped == 2)
        hom = pl.DataFrame({
            'id': var_ids[rows],
            'individual': sample_list[cols],
            'genotype': 2
        })
        geno_melt = pl.concat([geno_melt, hom])

    var_pheno_df = geno_melt.lazy().join(melted_pheno_df, on='individual', how='left')
    if debug:
        # Return intermediate dataframe if debugging
        return var_pheno_df
    
    var_pheno_df = var_pheno_df.group_by(
           ['id', 'phenotype']
           ).agg([
                pl.len().alias('n_individuals'),
                pl.col('pheno_value').mean().cast(pl.Float32).alias('mean_pheno_value'),
                pl.col('pheno_value').std().cast(pl.Float32).alias('std_pheno_value'),
            ]).drop_nulls(subset=['mean_pheno_value'])
    return var_pheno_df


anngeno_path = '/home/dnanexus/data_dir/genebass_1e6_coding_variants.ag'
eur_samples_path = '/home/dnanexus/data_dir/unrelated_cauc_samples_3rd_degree.csv'
# olink_path = '/home/dnanexus/data_dir/olink/protrider_lite_output/log2fc.csv'
olink_path = "/home/dnanexus/data_dir/olink/olink_corrected_rint_90_pcs.parquet"

sample_ids = zarr.open(f'{anngeno_path}/zarr_store/samples', mode='r')[:]
var_ids = pl.read_parquet(f'{anngeno_path}/variant_metadata.parquet', columns=['id'])['id'].to_numpy()
geno = zarr.open(f'{anngeno_path}/zarr_store/genotypes', mode='r')

eur_samples = pl.read_csv(eur_samples_path).rename({'eid': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
)['individual'].to_list()

olink_df = pl.read_parquet(olink_path).rename({'sample': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
).filter(
    pl.col("individual").is_in(eur_samples)
).fill_nan(None)

unique_phenotypes = [pheno for pheno in olink_df.columns if pheno != 'individual']

olink_melt = (
    olink_df.unpivot(
        index=['individual'],
        on=unique_phenotypes,
        variable_name='phenotype',
        value_name='pheno_value',
    )
    .with_columns(
        phenotype = pl.col('phenotype') + '_olink'
    )
    .lazy()
)

# Optimized approach to find indices of olink_samples in sample_ids
# Step A: Create a lookup dictionary mapping each ID in the large array to its original index. This takes O(N) time, where N is the size of sample_ids.
olink_sids = olink_df['individual'].to_numpy()
sample_id_to_index = {sid: i for i, sid in enumerate(sample_ids)}

# Step B: Iterate through the smaller array and find the index for each element if it exists in our lookup dictionary. This takes O(M) time, where M is the size of olink_sids.
found_indices = []
for sid in olink_sids:
    if sid in sample_id_to_index:
        found_indices.append(sample_id_to_index[sid])
print(f"Found {len(found_indices)} matching IDs.")

# Step C: Convert the list of indices to a NumPy array and sort it.
# Sorting ensures the output is identical to the original np.where approach, which returns indices in ascending order.
olink_indices = np.sort(found_indices)

output_dir = "/home/dnanexus/data_dir/olink_appv_chunks_EUR"
output_dir = "/home/dnanexus/data_dir/olink_appv_chunks_EUR_coding"
chunk_size = 10_000

for chunk_num in tqdm(range(var_ids.shape[0]//chunk_size + 1)):
    process_genotype_chunk(
        geno=geno[chunk_num*chunk_size:(chunk_num+1)*chunk_size, olink_indices],
        var_ids=var_ids[chunk_num*chunk_size:(chunk_num+1)*chunk_size],
        sample_list=sample_ids,
        melted_pheno_df=olink_melt,
        homozygous=False,
    ).sink_parquet(f"{output_dir}/variant_pheno_chunk{chunk_num}.parquet")

# Concat all files into one
files = [f"{output_dir}/variant_pheno_chunk{i}.parquet" for i in range(var_ids.shape[0]//chunk_size + 1)]
lazy_frames = [pl.scan_parquet(f) for f in files]
combined = pl.concat(lazy_frames)

combined.sink_parquet(f"{output_dir}/variant_pheno_EUR.parquet", engine='streaming')