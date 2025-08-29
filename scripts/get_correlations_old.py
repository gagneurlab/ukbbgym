import sys
import yaml
import zarr
import pandas as pd
import numpy as np
from tqdm import tqdm
from anngeno import AnnGeno
import statsmodels.api as sm
from joblib import Parallel, delayed


# Get covariate corrected phenotypes
def cov_prs_correction(all_df, phenotypes, covariates, prs_pheno_map):
    # Initialize an empty DataFrame to store residuals
    # all_df.set_index('sample', inplace=True)
    cov_prs_corrected_phenos = pd.DataFrame(
        index=all_df.index
    )  # Index is the sample ID

    # Perform linear regression for each phenotype
    for pheno in tqdm(phenotypes):
        # Drop NaN values for the current phenotype
        combined_df = all_df[[pheno] + covariates + [prs_pheno_map[pheno]]].dropna()
        y = combined_df[pheno]
        X = combined_df.drop(columns=[pheno])
        X = sm.add_constant(X)  # Add a constant term for the intercept

        # Fit the model
        model = sm.OLS(y, X).fit()

        # Save residuals
        residuals = pd.Series(model.resid, index=combined_df.index, name=pheno)
        cov_prs_corrected_phenos = pd.concat(
            [cov_prs_corrected_phenos, residuals], axis=1
        )

    # Reset the index for the resulting DataFrame
    cov_prs_corrected_phenos.reset_index(inplace=True)
    # cov_prs_corrected_phenos.columns = ['sample'] + [f"{pheno}_cov_prs_corrected" for pheno in phenotypes]
    return cov_prs_corrected_phenos


def pheno_burden_correlation(assoc_df, gt_df, annotation, correlation_type):
    """
    Computes gene-phenotype correlations for a given annotation.
    This function is parallelized using joblib for trait-gene pairs, focusing on associations.

    Args:
        assoc_df (pd.DataFrame): DataFrame with associations (phenotype, gene).
        gt_df (pd.DataFrame): DataFrame with gene burdens and corrected phenotypes.
        annotation (str): The current annotation being processed.
        correlation_type (str): Type of correlation to compute (e.g., 'spearman', 'pearson').

    Returns:
        pd.DataFrame: DataFrame containing gene-phenotype correlations.
    """

    # Define a helper function to process each trait-gene pair
    def _process_trait_gene(trait, gene, gt_df_local, corr_type, anno):
        """
        Calculates correlation for a single trait-gene pair.
        This function will be called in parallel.
        """
        pheno = trait.replace(" ", "_")  # Format phenotype name for column lookup
        correlation = np.nan
        # correlation_non_zero = np.nan # This was commented out in original, keeping it that way

        try:
            # Calculate correlation for the gene and phenotype
            correlation = (
                gt_df_local[[gene, pheno]].dropna().corr(method=corr_type).iloc[0, 1]
            )

        except ValueError:
            print(
                f"Wrong correlation type specified for {anno}, {pheno}, {gene}. Reverting to spearman."
            )
            correlation = (
                gt_df_local[[gene, pheno]].dropna().corr(method="spearman").iloc[0, 1]
            )
        except Exception as e:
            print(f"Cannot compute correlation for {anno}, {pheno}, {gene}. Error: {e}")
            correlation = np.nan

        # Return a DataFrame for the current trait-gene correlation
        return pd.DataFrame(
            {
                "annotation": anno,
                "phenotype": trait,
                "gene": gene,
                "correlation": correlation,
                # "correlation_non_zero": correlation_non_zero,
            },
            index=[0],
        )

    # Generate a list of all (trait, gene) pairs to process
    # This creates the iterable for joblib.Parallel
    tasks = []
    for trait in assoc_df.phenotype.unique():
        gene_list = list(assoc_df.query("phenotype == @trait")["region"].astype(str))
        for gene in gene_list:
            tasks.append((trait, gene))

    # Use joblib.Parallel to distribute the processing of each (trait, gene) pair. Pass a local copy of gt_df for each task to avoid potential shared memory issues if the original gt_df is modified (though not expected here).
    results = Parallel(n_jobs=-1)(
        delayed(_process_trait_gene)(trait, gene, gt_df, correlation_type, annotation)
        for trait, gene in tqdm(tasks, desc=f"Correlating for {annotation}")
    )

    # Concatenate all results into a single DataFrame
    return pd.concat(results)


def compute_correlations(
    config_path,
    zarr_burdens_path,
    genes_to_keep=None,
    max_burden=False,
    correlation_type="spearman",
):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    anngeno_file = config.get("anngeno_file")
    # annotation_cats = config.get("rare_variant_annotations")
    associations_df_path = config.get("associations_df_path")
    phenotypes = config.get("phenotypes_for_testing")
    covs = config.get("covariates")
    prs_pheno_map_file = config.get("prs_pheno_map_file")
    prs_file = config.get("prs_file")

    cov_pheno_df = pd.read_parquet(
        f"{anngeno_file}/phenotypes.parquet", columns=["sample"] + covs + phenotypes
    ).set_index("sample")
    prs_df = pd.read_parquet(prs_file)
    prs_df = prs_df[prs_df.index.isin(cov_pheno_df.index)]
    prs_pheno_map = pd.read_csv(prs_pheno_map_file)
    prs_pheno_map = dict(zip(prs_pheno_map["phenotype"], prs_pheno_map["pgs_id"]))
    all_df = pd.concat([cov_pheno_df, prs_df], axis=1)
    pheno_corrected_df = cov_prs_correction(all_df, phenotypes, covs, prs_pheno_map)

    zarr_group = zarr.open_group(zarr_burdens_path, mode="r")
    sample_list = zarr_group["samples"][:]
    gene_list = zarr_group["genes"][:]
    annotation_list = zarr_group["annotations"][:]

    assoc_df = pd.read_parquet(associations_df_path)
    if genes_to_keep is not None:
        assoc_df = assoc_df[assoc_df.gene.isin(genes_to_keep)]

    rho_df_sum_list = []
    rho_df_max_list = []
    print(f"Starting correlation computation for {len(annotation_list)} annotations")
    print(annotation_list)
    for anno in tqdm(annotation_list):
        anno_idx = np.where(annotation_list == anno)[0][0]
        sum_burdens = zarr_group["sum_burdens"][:, :, anno_idx]
        gt_df_sum = pd.DataFrame(
            sum_burdens, index=sample_list, columns=gene_list
        ).merge(pheno_corrected_df, left_index=True, right_on="sample")
        rho_df_sum_list.append(
            pheno_burden_correlation(assoc_df, gt_df_sum, anno, correlation_type)
        )
        if max_burden:
            max_burdens_zarr = zarr_group["max_burdens"][:, :, anno_idx]
            gt_df_max = pd.DataFrame(
                max_burdens_zarr, index=sample_list, columns=gene_list
            ).merge(pheno_corrected_df, left_index=True, right_on="sample")
            rho_df_max_list.append(
                pheno_burden_correlation(assoc_df, gt_df_max, anno, correlation_type)
            )

    rho_df_sum = pd.concat(rho_df_sum_list)
    rho_df_sum["aggregation"] = "sum"

    if max_burden:
        rho_df_max = pd.concat(rho_df_max_list)
        rho_df_max["aggregation"] = "max"
        rho_df = pd.concat([rho_df_sum, rho_df_max])
        return rho_df

    return rho_df_sum
