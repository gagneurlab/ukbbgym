import os
import sys
import yaml
import zarr
import pandas as pd
import polars as pl
import numpy as np
from tqdm import tqdm
from anngeno import AnnGeno

import click

def get_gene_burdens(
    anngeno_obj, gene_num, annotation_list, max_burden=False
):
    try:
        gene = anngeno_obj.get_region(gene_num)
    except Exception as e:
        print(f"Cannot find {gene_num} in AnnGeno, Returning NaNs. Error: {e}")
        nan_gis = np.zeros((len(anngeno_obj.samples), len(annotation_list))) * np.nan
        return nan_gis, nan_gis
    
    geno = gene["genotypes"]
    no_variant_mask = geno.sum(axis = 1) == 0

    try:
        var_scores = gene["annotations"][annotation_list].fillna(0).to_numpy().astype(np.float32)
    except Exception as e:
        print(f"Cannot find {annotation_list} in AnnGeno, Returning NaNs. Error: {e}")
        nan_gis = np.zeros((len(anngeno_obj.samples), len(annotation_list))) * np.nan
        return nan_gis, nan_gis
    
    # Calculate sum burden directly
    gis_sum = np.dot(geno, var_scores)
    gis_sum[no_variant_mask, :] = np.nan

    # If max_burden is False, return sum burden
    if not max_burden:
        return gis_sum, np.nan

    # For max calculation
    gis_max = []
    for a in tqdm(range(var_scores.shape[1])):
        gis_max.append(np.max(geno*var_scores[:, a], axis=1))
    gis_max = np.stack(gis_max, axis=1)
    gis_max[no_variant_mask, :] = np.nan

    return gis_sum, gis_max


def get_burdens_array(
    config,
    associations_df_path,
    annotation_list,
    new_anno_df=None,
    max_burden=False,
    only_snps=False,
    debug=False
):
    maf = config.get("association_testing_maf")

    associations_df = pl.read_parquet(associations_df_path)
    if debug:
        print("Debug is True, using only 5 associations")
        associations_df = associations_df.head()
    genes = associations_df['gene'].unique()

    print("Loading AnnGeno file")
    anngeno_file = config.get("anngeno_file")
    ag = AnnGeno(filename=anngeno_file, filemode="r")

    print(f"Filtering for variants with MAF < {maf}")

    # --- MODIFIED LINE START ---
    # Convert Polars Series to a Python list for direct use in isin() with pandas
    genes_list = genes.to_list()
    variants_to_keep_df = ag.annotations[
        (ag.annotations['MAF'] < maf) & (ag.annotations['region'].isin(genes_list))
    ]
    variants_to_keep = set(variants_to_keep_df["id"])
    # --- MODIFIED LINE END ---

    ag.subset_variants(variants_to_keep)

    if only_snps:
        print(f"Filtering for SNPs only")
        # Ensure that `ag.annotations` is still a pandas DataFrame after subsetting
        snp_variants = set(ag.annotations.query("(ref.str.len()==1) & (alt.str.len()==1)")["id"])
        ag.subset_variants(snp_variants)

    if new_anno_df is not None:
        # TODO improve this (subset anngeno)
        ## CAUTION very hacky and not stable
        print("Adding annotations to AnnGeno file")
        # new_cols = list(set(annos.columns) - set(ag.annotations.columns))
        new_cols = list(set(new_anno_df.columns) - set(ag.annotations.columns))
        print(f"Adding new annotations to anngeno: {new_cols}")
        merged = ag.annotations\
            .merge(new_anno_df[["chrom", "pos", "ref", "alt", "region", *new_cols]], how = "left", on = ["chrom", "pos", "ref", "alt", "region"])
        ag._set_annotations(merged)

    gene_id_list = list(genes)
    gene_burdens_sum = []
    gene_burdens_max = []
    print(f"Starting to compute gene burdens for {len(genes)} genes")
    for gene in tqdm(gene_id_list):
        gis_sum, gis_max = get_gene_burdens(ag, gene, annotation_list, max_burden)
        gene_burdens_sum.append(gis_sum)
        if max_burden:
            gene_burdens_max.append(gis_max)

    gene_burdens_sum_df = np.stack(gene_burdens_sum, axis=1)  # (n_samples, n_genes, n_annotations)
    if max_burden:
        gene_burdens_max_df = np.stack(gene_burdens_max, axis=1)  # (n_samples, n_genes, n_annotations)
        print("Returning sum and max burden.")
        return gene_burdens_sum_df, gene_burdens_max_df, ag.samples, gene_id_list

    print("Max burden is False, returning only sum burden.")
    return gene_burdens_sum_df, gene_burdens_sum_df, ag.samples, gene_id_list


# Example usage
# To compute burdens based on the config's rare_variant_annotations:
# compute_and_store_burdens(
#     config_path="your_config.yaml",
#     output_zarr="output.zarr",
#     max_burden=True,
#     only_snps=True,
#     overwrite=False,
# )

# To add new annotations from an annotation scores file:
# compute_and_store_burdens(
#     config_path="your_config.yaml",
#     output_zarr="output.zarr",
#     anno_scores_path="annotation_scores.parquet",
#     max_burden=False,
#     overwrite=False,
# )

def compute_and_store_burdens(
    config_path,
    output_zarr,
    anno_scores_path=None,
    max_burden=False,
    only_snps=False,
    overwrite=False,
    debug=False,
):
    """
    Computes and stores variant burdens in a Zarr array, handling both initial creation
    and adding new annotations.

    Args:
        config_path (str): Path to the configuration YAML file.
        output_zarr (str): Path to the output Zarr file.
        anno_scores_path (str, optional): Path to a Parquet file containing annotation scores.
                                          If provided, new annotations from this file will be added.
                                          Defaults to None.
        max_burden (bool, optional): Whether to compute and store the maximum burden.
                                     Defaults to False.
        only_snps (bool, optional): Whether to consider only SNPs for burden calculation.
                                   Defaults to False.
        overwrite (bool, optional): Whether to overwrite the existing Zarr file.
                                    Defaults to False.
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    associations_df_path = config.get("associations_df_path")
    all_annotation_list = []

    if anno_scores_path:
        anno_scores_df = pl.read_parquet(anno_scores_path)
        available_annotations = list(set(anno_scores_df.columns) - set(["chrom", "pos", "ref", "alt", "region"]))
        all_annotation_list.extend(available_annotations)
    else:
        rare_variant_annotations_dict = config.get('rare_variant_annotations')
        if rare_variant_annotations_dict:
            for category in rare_variant_annotations_dict.values():
                all_annotation_list.extend(category)

    zarr_file_path = output_zarr
    anno_chunk_size = 1

    if os.path.exists(zarr_file_path):
        if overwrite:
            print(f"Overwriting existing Zarr file at {zarr_file_path}.")
            try:
                zarr.rmtree(zarr_file_path)
            except Exception as e:
                print(f"Error deleting existing Zarr file: {e}")
                sys.exit(1)
        else:
            print(f"Zarr file exists at {zarr_file_path}, checking for new annotations.")
            try:
                root = zarr.group(zarr_file_path, mode='r+')
                sum_burdens = root.get("sum_burdens")
                zarr_annotations = root.get("annotations")
                max_burdens = root.get("max_burdens")

                if sum_burdens is None or zarr_annotations is None:
                    print("Existing Zarr file is incomplete. Consider overwriting or creating a new one.")
                    sys.exit(1)

                existing_annotations = list(zarr_annotations[:])
                new_annotation_list = list(set(all_annotation_list) - set(existing_annotations))
                union_annotation_list = existing_annotations + new_annotation_list

                if new_annotation_list:
                    print(f"New annotations found: {new_annotation_list}")
                    get_burdens_kwargs = {
                        "config": config,
                        "associations_df_path": associations_df_path,
                        "annotation_list": new_annotation_list,
                        "max_burden": max_burden,
                        "only_snps": only_snps,
                        "debug": debug,
                    }
                    if anno_scores_path:
                        get_burdens_kwargs["new_anno_df"] = anno_scores_df

                    new_gene_burdens_sum_df, new_gene_burdens_max_df, _, _ = get_burdens_array(**get_burdens_kwargs)

                    current_shape = sum_burdens.shape
                    new_shape = (current_shape[0], current_shape[1], current_shape[2] + len(new_annotation_list))
                    sum_burdens.resize(new_shape)
                    sum_burdens[:, :, current_shape[2]:] = new_gene_burdens_sum_df

                    if max_burden and max_burdens is not None:
                        print("Max burden is True, appending max_burdens array.")
                        max_burdens.resize(new_shape)
                        max_burdens[:, :, current_shape[2]:] = new_gene_burdens_max_df

                    if zarr_annotations is not None:
                        zarr_annotations_new_shape = (len(union_annotation_list),)
                        if zarr_annotations.shape != zarr_annotations_new_shape:
                            zarr_annotations.resize(zarr_annotations_new_shape)
                        zarr_annotations[:] = union_annotation_list
                        print("Updated 'annotations' array.")
                    print(f"Appended data for new annotations: {new_annotation_list}")
                else:
                    print("No new annotations to add.\nExiting.")

            except Exception as e:
                print(f"Error accessing or updating existing Zarr file: {e}")
                sys.exit(1)

    else:
        print(f"Zarr file does not exist at {zarr_file_path}, creating a new one.")
        get_burdens_kwargs = {
            "config": config,
            "associations_df_path": associations_df_path,
            "annotation_list": all_annotation_list,
            "max_burden": max_burden,
            "only_snps": only_snps,
            "debug": debug,
        }
        if anno_scores_path:
            get_burdens_kwargs["new_anno_df"] = anno_scores_df

        gene_burdens_sum_df, gene_burdens_max_df, sample_id_arr, gene_id_list = get_burdens_array(**get_burdens_kwargs)

        try:
            root = zarr.group(zarr_file_path)
            root.create_array(
                "sum_burdens",
                shape=gene_burdens_sum_df.shape,
                dtype=gene_burdens_sum_df.dtype,
                chunks=(gene_burdens_sum_df.shape[0], gene_burdens_sum_df.shape[1], anno_chunk_size),
                overwrite=overwrite,
            )[:] = gene_burdens_sum_df

            if max_burden:
                root.create_array(
                    "max_burdens",
                    shape=gene_burdens_max_df.shape,
                    dtype=gene_burdens_max_df.dtype,
                    chunks=(gene_burdens_max_df.shape[0], gene_burdens_max_df.shape[1], anno_chunk_size),
                    overwrite=overwrite,
                )[:] = gene_burdens_max_df

            root.create_array("samples", shape=sample_id_arr.shape, dtype="str", overwrite=overwrite)[:] = sample_id_arr
            root.create_array("genes", shape=len(gene_id_list), dtype="str", overwrite=overwrite)[:] = gene_id_list
            root.create_array("annotations", shape=len(all_annotation_list), dtype="str", overwrite=overwrite)[:] = all_annotation_list
            print("Created new zarr array 'sum_burdens', 'samples', 'genes', and 'annotations'.")

        except Exception as e:
            print(f"Error creating new Zarr file: {e}")
            sys.exit(1)


@click.group()
def cli():
    pass

@cli.command()
@click.option("--config-path", type=str, required=True, help="Config file with all details")
@click.option("--output-zarr", type=str, required=True, help="Output zarr file")
@click.option("--max-burden", is_flag=True, default=False, help="Compute burdens using max scores")
@click.option("--only-snps", is_flag=True, default=False, help="Use only SNPs to compute burdens")
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite old zarr file") #TODO
@click.option("--debug", is_flag=True, default=False, help="Use only 5 associaitons to debug code")
def compute_burdens(
    config_path: str,
    output_zarr: str,
    max_burden: bool = False,
    only_snps: bool = False,
    overwrite: bool = False,
    debug: bool = False,
):
    print('You are running the script to compute gene burdens')

    compute_and_store_burdens(
        config_path=config_path,
        output_zarr=output_zarr,
        max_burden=max_burden,
        only_snps=only_snps,
        overwrite=overwrite,
        debug=debug
    )

    print('Gene burdens have been computed and stored')

if __name__ == "__main__":
    cli()