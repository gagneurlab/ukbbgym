import os
import yaml
import pandas as pd
import polars as pl
import statsmodels.api as sm
from tqdm import tqdm


def cov_prs_correction(all_df, phenotypes, covariates=None, prs_pheno_map=None):
    cov_prs_corrected_phenos = pd.DataFrame(index=all_df.index)

    for pheno in tqdm(phenotypes, desc="Correcting phenotypes"):
        if (prs_pheno_map is not None) and (covariates is not None):
            combined_df = all_df[[pheno] + covariates + [prs_pheno_map[pheno]]].dropna()
        elif covariates is not None:
            combined_df = all_df[[pheno] + covariates].dropna()
        elif prs_pheno_map is not None:
            combined_df = all_df[[pheno] + [prs_pheno_map[pheno]]].dropna()
        else:
            return all_df[[pheno]].dropna()

        y = combined_df[pheno]
        X = combined_df.drop(columns=[pheno])
        X = sm.add_constant(X)

        model = sm.OLS(y, X).fit()
        residuals = pd.Series(model.resid, index=combined_df.index, name=pheno)
        cov_prs_corrected_phenos = pd.concat([cov_prs_corrected_phenos, residuals], axis=1)

    cov_prs_corrected_phenos.reset_index(inplace=True)
    return cov_prs_corrected_phenos


def compute_correlations_lazy(df_lazy: pl.LazyFrame) -> pl.LazyFrame:
    df_clean = df_lazy.filter(
        pl.col("value").is_not_nan() &
        pl.col("value").is_not_null() &
        pl.col("value").is_finite()
    )

    df_ranks = df_clean.with_columns([
        pl.col("sum").rank().over(["annotation", "phenotype"]).alias("rank_sum"),
        pl.col("max").rank().over(["annotation", "phenotype"]).alias("rank_max"),
        pl.col("top2").rank().over(["annotation", "phenotype"]).alias("rank_top2"),
        pl.col("value").rank().over(["annotation", "phenotype"]).alias("rank_value"),
    ])

    correlations = df_ranks.group_by(["annotation", "phenotype"]).agg([
        pl.corr("rank_sum", "rank_value").alias("sum_spearman"),
        pl.corr("rank_max", "rank_value").alias("max_spearman"),
        pl.corr("rank_top2", "rank_value").alias("top2_spearman"),
    ])

    del df_ranks
    return correlations


def compute_correlations(
    burdens_dir: str,
    pheno_corrected_df: pl.DataFrame,
    assocs_df: pl.DataFrame,
    config: dict,
    filter_nan: bool = True,
    subset_annos: list | None = None,
    subset_samples: list | None = None,
) -> pl.DataFrame:

    pheno_df = pheno_corrected_df.unpivot(
        index=['sample_id'],
        variable_name='phenotype',
        value_name='value'
    )

    corr_df_list = []

    for gene_id in tqdm(assocs_df['gene_id'].unique(), desc="Correlations for genes"):
        try:
            bdf = pl.scan_parquet(f'{burdens_dir}/{gene_id}.parquet')
        except FileNotFoundError:
            print(f"File for gene {gene_id} not found. Skipping.")
            continue

        if subset_annos:
            bdf = bdf.filter(pl.col('annotation').is_in(subset_annos))
        if subset_samples:
            bdf = bdf.filter(pl.col('sample_id').is_in(subset_samples))

        if filter_nan:
            bdf_filtered = bdf.filter(
                pl.all_horizontal(pl.col(['sum', 'max', 'top2']).is_not_nan())
            )
        else:
            cols_to_fill = ['max', 'sum', 'top2']
            modes = bdf.group_by("annotation").agg([
                pl.col(col).drop_nans().mode().first().alias(f"{col}_mode")
                for col in cols_to_fill
            ])
            bdf_with_modes = bdf.join(modes, on="annotation")
            bdf_filtered = bdf_with_modes.with_columns([
                pl.when(pl.col(col).is_nan())
                .then(pl.col(f"{col}_mode"))
                .otherwise(pl.col(col))
                .alias(col)
                for col in cols_to_fill
            ]).drop([f"{col}_mode" for col in cols_to_fill])

        phenos_needed = assocs_df.filter(pl.col('gene_id') == gene_id)['phenotype'].unique().to_list()
        pheno_filtered = pheno_df.filter(pl.col('phenotype').is_in(phenos_needed)).lazy()
        cdf = bdf_filtered.join(pheno_filtered, on='sample_id', how='left')

        corr_df_list.append(
            compute_correlations_lazy(cdf).with_columns(
                pl.lit(gene_id).alias('gene_id')
            ).collect()
        )

    corr_df = pl.concat(corr_df_list)

    rare_variant_annotations_dict = config.get('rare_variant_annotations', {})
    annotation_category_map = {
        ann: category
        for category, anns in rare_variant_annotations_dict.items() if category != 'misc'
        for ann in anns
    }

    corr_df = corr_df.with_columns(
        pl.col("annotation").replace(annotation_category_map).alias("category")
    )

    corr_long = corr_df.unpivot(
        index=["annotation", "phenotype", "gene_id", "category"],
        variable_name="correlation_type",
        value_name="correlation"
    ).with_columns([
        pl.col("correlation_type").str.extract(r"(pearson|spearman)").alias("method"),
        pl.col("correlation_type").str.extract(r"(sum|max|top2)").alias("aggregation"),
    ])

    del corr_df
    return corr_long


# def compute_correlations_wrapper(
#     config_path: str,
#     burdens_dir: str,
#     associations_file: str,
#     pheno_file: str,
#     subset_samples: list | None = None,
#     subset_annos: list | None = None,
#     filter_nan: bool = False,
# ) -> pl.DataFrame:

#     with open(config_path) as f:
#         config = yaml.safe_load(f)

#     covs = config.get("covariates")

#     assocs_df = pl.read_parquet(associations_file)
#     phenotypes = list(assocs_df['phenotype'].unique())
#     pheno_list = [p + '_prs_corrected' for p in phenotypes]

#     pheno_df = pl.read_parquet(pheno_file, columns=["eid"] + covs + pheno_list).rename({'eid': 'sample_id'})
#     corrected_df = cov_prs_correction(
#         pheno_df.to_pandas().set_index('sample_id'),
#         pheno_list,
#         covariates=covs
#     )
#     pheno_corrected_df = pl.from_pandas(corrected_df)
#     pheno_corrected_df.columns = ['sample_id'] + phenotypes

#     corr_long = compute_correlations(
#         burdens_dir=burdens_dir,
#         pheno_corrected_df=pheno_corrected_df,
#         assocs_df=assocs_df,
#         config=config,
#         filter_nan=filter_nan,
#         subset_annos=subset_annos,
#         subset_samples=subset_samples,
#     )

#     return corr_long
