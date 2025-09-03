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

from anngeno import AnnGeno
import multiprocessing

import matplotlib.pyplot as plt
from plotnine import *

num_cores = multiprocessing.cpu_count()
print(num_cores)
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


pheno_path = '/home/dnanexus/data_dir/phenotypes/phenotypes190_missing20_unique2_int.parquet'
eur_samples_path = '/home/dnanexus/data_dir/unrelated_cauc_samples_3rd_degree.csv'
anngeno_path = '/home/dnanexus/data_dir/genebass_1e6_coding_variants.ag'
output_dir = "/home/dnanexus/data_dir/var_pheno_chunks"

sample_ids = zarr.open(f'{anngeno_path}/zarr_store/samples', mode='r')[:]
var_ids = pl.read_parquet(f'{anngeno_path}/variant_metadata.parquet', columns=['id'])['id'].to_numpy()

geno = zarr.open(f'{anngeno_path}/zarr_store/genotypes', mode='r')

eur_samples = pl.read_csv(eur_samples_path).rename({'eid': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
)['individual'].to_list()

pheno = pl.read_parquet(pheno_path).drop(['FID']).rename({'IID': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
).filter(
    pl.col("individual").is_in(eur_samples)
).unpivot(
    index=['individual'],
    variable_name='phenotype',
    value_name='pheno_value',
).drop_nans().drop_nulls().lazy()

chunk_size = 10_000
for chunk_num in tqdm(range(var_ids.shape[0]//chunk_size + 1)):
    process_genotype_chunk(
        geno=geno[chunk_num*chunk_size:(chunk_num+1)*chunk_size],
        var_ids=var_ids[chunk_num*chunk_size:(chunk_num+1)*chunk_size],
        sample_list=sample_ids,
        melted_pheno_df=pheno,
        homozygous=False,
    ).sink_parquet(f"{output_dir}/variant_pheno_chunk{chunk_num}.parquet")