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

num_cores = multiprocessing.cpu_count()
print(num_cores)
def process_phenotypes_prs_long(
    pheno_lazy: str,
    prs_lazy: str,
    cov_lazy: str,
    unique_phenotypes: pl.DataFrame,
    cov_list: list,
    quantitative: bool = True
) -> pl.DataFrame:
    """
    Process phenotypes and PRS, compute residuals for each phenotype, and return a long-format Polars DataFrame.
    """
    print("Process phenotypes and PRS, compute residuals for each phenotype, and return a long-format Polars DataFrame.")

    # --- Merge all into one lazy DataFrame ---
    all_lazy = pheno_lazy.join(prs_lazy, on='individual', how='inner').join(cov_lazy, on='individual', how='inner')
    
    # Collect once for regression computations (still needed for statsmodels)
    all_pd = all_lazy.collect().to_pandas()

    # --- Compute residuals ---
    all_residuals_dfs = []

    for phenotype in tqdm(unique_phenotypes):
        pheno_cols = [phenotype, f"{phenotype}_prs"] + cov_list
        temp_df = all_pd[['individual'] + pheno_cols].dropna()
        if len(temp_df) == 0:
            print(f"No data for phenotype: {phenotype}")
            continue

        y = temp_df[phenotype]
        X = temp_df.drop(columns=[phenotype, 'individual'])
        X = sm.add_constant(X)

        if quantitative:
            model = sm.OLS(y, X).fit()
            residuals = pd.Series(model.resid, index=temp_df.index, name=f"{phenotype}_residual")
        else:
            model = sm.GLM(y, X, family=sm.families.Binomial()).fit()
            residuals = pd.Series(model.resid_deviance, index=temp_df.index, name=f"{phenotype}_residual")

        pheno_residuals = pd.concat([temp_df[['individual']], residuals], axis=1)
        all_residuals_dfs.append(pheno_residuals)


    if not all_residuals_dfs:
        raise ValueError("No residuals could be computed")

    # --- Step 5: Convert to long-format lazy DataFrame ---
    long_lazy_dfs = []
    for residual_df in all_residuals_dfs:
        p_wide_lazy = pl.LazyFrame(residual_df).with_columns(
            pl.col('individual').cast(pl.String)
        )
        pheno_cols = [c for c in residual_df.columns if c.endswith('_residual')]
        pdf_lazy = (
            p_wide_lazy.unpivot(
                index=['individual'],
                on=pheno_cols,
                variable_name='phenotype',
                value_name='pheno_value',
            )
            .with_columns(
                pl.col('phenotype').str.replace('_residual', '').alias('phenotype')
            )
        )
        long_lazy_dfs.append(pdf_lazy)

    combined_pdf_lazy = pl.concat(long_lazy_dfs) if len(long_lazy_dfs) > 1 else long_lazy_dfs[0]

    return combined_pdf_lazy
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

pheno_path = 'PATH_TO_FILE'
prs_path = 'PATH_TO_FILE'
cov_path = 'PATH_TO_FILE'

anngeno_path = 'PATH_TO_FILE'
eur_samples_path = 'PATH_TO_FILE'
sample_ids = zarr.open(f'{anngeno_path}/zarr_store/samples', mode='r')[:]
sample_ids
var_ids = pl.read_parquet(f'{anngeno_path}/variant_metadata.parquet', columns=['id'])['id'].to_numpy()
var_ids
geno = zarr.open(f'{anngeno_path}/zarr_store/genotypes', mode='r')
geno[:5].shape
eur_samples = pl.read_csv(eur_samples_path).rename({'eid': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
)['individual'].to_list()

pheno_lazy = pl.read_parquet(pheno_path).drop(['FID']).rename({'IID': 'individual'}).with_columns(
    pl.col("individual").cast(pl.Utf8)
).filter(
    pl.col("individual").is_in(eur_samples)
).fill_nan(None).lazy()

unique_phenotypes = [pheno for pheno in pheno_lazy.columns if pheno != 'individual']
quant_phenotypes = [pheno for pheno in unique_phenotypes if "jurgens" not in pheno]

prs_cols = [f"{pheno}_prs" for pheno in quant_phenotypes]
prs_lazy = (
    pl.scan_parquet(prs_path)
    .fill_nan(None)
    .rename({'IID': 'individual'})
    .select(['individual'] + prs_cols)
    .with_columns(pl.col("individual").cast(pl.Utf8))
    .drop_nulls()
)

config_path = 'PATH_TO_FILE'
with open(config_path) as f:
    config = yaml.safe_load(f)

cov_list = config.get("covariates", [])
cov_lazy = (
    pl.scan_parquet(cov_path)
    .rename({'sample': 'individual'})
    .select(['individual'] + cov_list)
)

corr_phenos_lazy = process_phenotypes_prs_long(
    pheno_lazy=pheno_lazy,
    prs_lazy=prs_lazy,
    cov_lazy=cov_lazy,
    unique_phenotypes=quant_phenotypes,
    cov_list=cov_list,
    quantitative=True
)

corr_phenos_lazy.head().collect()
output_dir = "PATH_TO_FILE"
chunk_size = 10_000

for chunk_num in tqdm(range(var_ids.shape[0]//chunk_size + 1)):
    process_genotype_chunk(
        geno=geno[chunk_num*chunk_size:(chunk_num+1)*chunk_size],
        var_ids=var_ids[chunk_num*chunk_size:(chunk_num+1)*chunk_size],
        sample_list=sample_ids,
        melted_pheno_df=corr_phenos_lazy,
        homozygous=False,
    ).sink_parquet(f"{output_dir}/variant_pheno_chunk{chunk_num}.parquet")