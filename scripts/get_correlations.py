import sys
import yaml
import zarr
import numpy as np
import polars as pl
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression


# -----------------------------
# Covariate + PRS correction
# -----------------------------
def cov_prs_correction(all_df: pl.DataFrame, phenotypes, covariates, prs_pheno_map):
    corrected_cols = []

    for pheno in tqdm(phenotypes, desc="Covariate correction"):
        try:
            pheno_prs = prs_pheno_map[pheno]
            cols = [pheno] + covariates + [pheno_prs] + ["sample"]

            df = all_df.select(cols).drop_nulls()

            y = df[pheno]
            X = df.drop([pheno, "sample"]).to_numpy()
            reg = LinearRegression().fit(X, y)
            
            corrected_cols.append(pl.DataFrame({"sample": df["sample"], pheno: y - reg.predict(X)}))
        except Exception as e:
            print(f"Failed to correct phenotype '{pheno}': {e}")

    # Join all corrected phenotypes
    if corrected_cols:
        corrected_df = pl.concat(corrected_cols, how='align_full')
        return corrected_df
    else:
        return pl.DataFrame([])


# -----------------------------
# Correlation computation
# -----------------------------
def gene_pheno_correlation_lazy(assoc_df: pl.DataFrame, gt_df: pl.DataFrame, annotation: str, correlation_type: str = "spearman"):
    """
    Compute correlations between gene burdens and phenotypes using Polars LazyFrames.
    """
    results = []

    for pheno in tqdm(assoc_df["phenotype"].unique().to_list(), desc=f"Processing phenotypes for {annotation}"):
        pheno_col = pheno.replace(" ", "_")  # match corrected phenotype column
        genes = assoc_df.filter(pl.col("phenotype") == pheno)["region"].cast(str).unique().to_list()

        for gene in genes:
            try:
                # LazyFrame to compute correlation
                lazy_df = (
                    gt_df.lazy()
                    .select([gene, pheno_col])
                    .drop_nans()
                    .select([
                        pl.corr(gene, pheno_col, method=correlation_type).alias("correlation")
                    ])
                )
                corr = lazy_df.collect().item()  # get scalar
                results.append({"annotation": annotation, "phenotype": pheno, "gene": gene, "correlation": corr})
            except Exception as e:
                results.append({"annotation": annotation, "phenotype": pheno, "gene": gene, "correlation": np.nan})
                print(f"Failed correlation for {pheno} x {gene}: {e}")

    return pl.DataFrame(results)


def compute_correlations(
    config_path, 
    zarr_burdens_path,
    genes_to_keep=None,
    max_burden=False,
    correlation_type='spearman',
):
    correlation_type = correlation_type.lower()
    if correlation_type not in ["pearson", "spearman"]:
        raise ValueError(f"Invalid correlation type: {correlation_type}. Defaulting to 'spearman'.")
    
    with open(config_path) as f:
        config = yaml.safe_load(f)

    anngeno_file = config["anngeno_file"]
    phenotypes = config["phenotypes_for_testing"]
    covs = config["covariates"]
    prs_file = config["prs_file"]
    prs_pheno_map_file = config["prs_pheno_map_file"]
    associations_df_path = config["associations_df_path"]

    # Load data
    cov_pheno_df = pl.read_parquet(f"{anngeno_file}/phenotypes.parquet").select(["sample"] + covs + phenotypes)
    prs_df = pl.read_parquet(prs_file).filter(pl.col("sample").is_in(cov_pheno_df["sample"]))
    prs_pheno_map = pl.read_csv(prs_pheno_map_file).to_dict(as_series=False)
    prs_pheno_map = dict(zip(prs_pheno_map["phenotype"], prs_pheno_map["pgs_id"]))

    # Merge covariates and PRS
    all_df = cov_pheno_df.join(prs_df, on="sample", how="inner")

    # Covariate and PRS correction
    pheno_corrected_df = cov_prs_correction(all_df, phenotypes, covs, prs_pheno_map)

    # Load zarr burden data
    zarr_group = zarr.open_group(zarr_burdens_path, mode="r")
    sample_list = zarr_group["samples"][:]
    gene_list = zarr_group["genes"][:]
    annotation_list = zarr_group["annotations"][:]

    # Load associations
    assoc_df = pl.read_parquet(associations_df_path)
    if genes_to_keep is not None:
        assoc_df = assoc_df.filter(pl.col("region").is_in(genes_to_keep))

    rho_df_sum_list = []
    rho_df_max_list = []
    rho_df_top2_list = []
    print(f"Starting {correlation_type} correlation computation for {len(annotation_list)} annotations")
    print(annotation_list)
    for anno in tqdm(annotation_list):
        anno_idx = np.where(annotation_list == anno)[0][0]
        
        sum_burdens_zarr = zarr_group["sum_burdens"][:, :, anno_idx]
        sum_burden_df = pl.DataFrame(sum_burdens_zarr, schema=list(gene_list)).with_columns([
            pl.Series(name="sample", values=sample_list)
        ])
        gt_df_sum = sum_burden_df.join(pheno_corrected_df, on="sample", how="inner")
        rho_df_sum_list.append(gene_pheno_correlation_lazy(assoc_df, gt_df_sum, anno, correlation_type))

        if max_burden:
            max_burdens_zarr = zarr_group["max_burdens"][:, :, anno_idx]
            max_burdens_df = pl.DataFrame(max_burdens_zarr, schema=list(gene_list)).with_columns([
                pl.Series(name="sample", values=sample_list)
            ])
            gt_df_max = max_burdens_df.join(pheno_corrected_df, on="sample", how="inner")
            rho_df_max_list.append(gene_pheno_correlation_lazy(assoc_df, gt_df_max, anno, correlation_type))
            
            top2_burdens_zarr = zarr_group["top2_burdens"][:, :, anno_idx]
            top2_burdens_df = pl.DataFrame(top2_burdens_zarr, schema=list(gene_list)).with_columns([
                pl.Series(name="sample", values=sample_list)
            ])
            gt_df_top2 = top2_burdens_df.join(pheno_corrected_df, on="sample", how="inner")
            rho_df_top2_list.append(gene_pheno_correlation_lazy(assoc_df, gt_df_top2, anno, correlation_type))


    rho_df_sum = pl.concat(rho_df_sum_list)
    rho_df_sum = rho_df_sum.with_columns(pl.lit("sum").alias("aggregation"))

    if max_burden:
        rho_df_max = pl.concat(rho_df_max_list)
        rho_df_max = rho_df_max.with_columns(pl.lit("max").alias("aggregation"))

        rho_df_top2 = pl.concat(rho_df_top2_list)
        rho_df_top2 = rho_df_top2.with_columns(pl.lit("top2").alias("aggregation"))

        rho_df = pl.concat([rho_df_sum, rho_df_max, rho_df_top2])
        return rho_df

    return rho_df_sum
