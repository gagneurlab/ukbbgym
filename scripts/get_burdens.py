import os
import sys
import yaml
import zarr
import click
import shutil
import pandas as pd
import polars as pl
import numpy as np

import torch
from tqdm import tqdm
from anngeno import AnnGeno
from joblib import Parallel, delayed

from numba import njit, prange
import multiprocessing

def get_gene_burdens(
    region_genotypes,
    region_annotations,
    annotation_list, 
    max_burden=False
):

    no_variant_mask = region_genotypes.sum(axis = 0) == 0

    try:
        var_scores = region_annotations[annotation_list].fill_nan(0).to_numpy().astype(np.float32).transpose()  # shape: (annotations, variants)
    except Exception as e:
        print(f"Error: {e}\nReturning NaNs.")
        return np.nan, np.nan, np.nan

    # Calculate sum burden directly
    gis_sum = np.dot(var_scores, region_genotypes).transpose()  # shape: (samples, annotations)
    gis_sum[no_variant_mask, :] = np.nan

    # If max_burden is False, return sum burden
    if not max_burden:
        return gis_sum, np.nan, np.nan # Still return a tuple to maintain consistent return type
    
    gis_max_list = []
    gis_top2_sum_list = []
    for a in range(var_scores.shape[0]):
        burden = np.abs(np.expand_dims(var_scores[a, :], axis=1) * region_genotypes)  # shape: (variants, samples)

        # Get top-k values per sample
        top2 = np.partition(burden, -2, axis=0)[-2:, :]  # shape: (k, samples)

        # Compute max (top-1) and sum of top-k
        max_vals = np.max(top2, axis=0)
        top2_sum = np.sum(top2, axis=0)

        gis_max_list.append(max_vals)
        gis_top2_sum_list.append(top2_sum)

    gis_max = np.stack(gis_max_list, axis=0).transpose()    # shape: (samples, annotations)
    gis_top2 = np.stack(gis_top2_sum_list, axis=0).transpose()
    
    # Handle no-variant case
    gis_max[no_variant_mask, :] = np.nan
    gis_top2[no_variant_mask, :] = np.nan

    return gis_sum, gis_max, gis_top2


@njit(parallel=True)
def compute_max_and_top2_chunked(score_vec, region_genotypes, chunk_size):
    n_variants, n_samples = region_genotypes.shape
    n_chunks = (n_samples + chunk_size - 1) // chunk_size

    max_vals = np.empty(n_samples, dtype=np.float32)
    top2_sums = np.empty(n_samples, dtype=np.float32)

    for c in prange(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, n_samples)
        for s in range(start, end):
            burden = np.abs(score_vec * region_genotypes[:, s])
            if len(burden) >= 2:
                top2 = np.partition(burden, -2)[-2:]
                max_vals[s] = top2.max()
                top2_sums[s] = top2.sum()
            elif len(burden) == 1:
                max_vals[s] = burden[0]
                top2_sums[s] = burden[0]
            else:
                max_vals[s] = 0.0
                top2_sums[s] = 0.0

    return max_vals, top2_sums

def get_gene_burdens_numba(
    region_genotypes,
    region_annotations,
    annotation_list, 
    max_burden=False,
    chunk_size=None,
):
    no_variant_mask = region_genotypes.sum(axis = 0) == 0

    try:
        var_scores = region_annotations[annotation_list].fill_nan(0).to_numpy().astype(np.float32).transpose()  # shape: (annotations, variants)
    except Exception as e:
        print(f"Error: {e}\nReturning NaNs.")
        return np.nan, np.nan, np.nan

    # Calculate sum burden directly
    gis_sum = np.dot(var_scores, region_genotypes).transpose()  # shape: (samples, annotations)
    gis_sum[no_variant_mask, :] = np.nan
    
    # If max_burden is False, return sum burden
    if not max_burden:
        return gis_sum, np.nan, np.nan # Still return a tuple to maintain consistent return type
    
    # Determine chunk size if not provided
    if chunk_size is None:
        num_cores = multiprocessing.cpu_count()
        chunk_size = max(1, region_genotypes.shape[1] // (num_cores * 2))
    
    print(f"Numba: Computing max and top2sum using chunk size: {chunk_size}")
    # Compute max + top2 via numba
    gis_max_list = []
    gis_top2_sum_list = []
    for a in tqdm(range(var_scores.shape[0])):
        score_vec = var_scores[a, :]
        max_vals, top2_sum = compute_max_and_top2_chunked(score_vec, region_genotypes, chunk_size)
        gis_max_list.append(max_vals)
        gis_top2_sum_list.append(top2_sum)

    gis_max = np.stack(gis_max_list, axis=0).T         # (samples, annotations)
    gis_top2 = np.stack(gis_top2_sum_list, axis=0).T   # (samples, annotations)

    gis_max[no_variant_mask, :] = np.nan
    gis_top2[no_variant_mask, :] = np.nan

    return gis_sum, gis_max, gis_top2

def get_gene_burdens_torch(
    region_genotypes,
    region_annotations,
    annotation_list,
    max_burden=False,
    device="cuda"
):
    with torch.no_grad():
        G = torch.tensor(region_genotypes, dtype=torch.float32, device=device)  # (variants, samples)
        A = torch.tensor(region_annotations[annotation_list].fill_nan(0).to_numpy(), dtype=torch.float32, device=device).transpose(0,1)  # (annotations, variants)

        no_variant_mask = G.sum(dim=0) == 0
        gis_sum = (A @ G).transpose(0,1) # shape: (samples, annotations)
        gis_sum[no_variant_mask] = float('nan')

        if not max_burden:
            return gis_sum.cpu().numpy(), None, None

        gis_max = []
        gis_top2 = []

        for a in range(A.shape[0]):
            scores = torch.abs(A[a, :] * G)  # (variants, samples)
            scores[G == 0] = float('-inf')

            top2_vals, _ = torch.topk(scores, k=2, dim=0, largest=True, sorted=False)  # (k, samples)
            gis_max.append(torch.max(top2_vals, dim=0).values)
            gis_top2.append(top2_vals.sum(dim=0))

            torch.cuda.empty_cache()
            
        gis_max = torch.stack(gis_max, dim=0).transpose(0,1)   # (samples, annotations)
        gis_top2 = torch.stack(gis_top2, dim=0).transpose(0,1) # (samples, annotations)

        gis_max[no_variant_mask] = float('nan')
        gis_top2[no_variant_mask] = float('nan')

        return gis_sum.cpu().numpy(), gis_max.cpu().numpy(), gis_top2.cpu().numpy()

def get_burdens_array(
    anngeno_path,   
    associations_df_path,
    maf,
    annotation_list,
    new_anno_df=None,
    max_burden=False,
    only_snps=False,
    debug=False,
    batch_size=32,
    n_jobs=32,
    device="cuda" if torch.cuda.is_available() else "cpu",
):
    print("Loading AnnGeno file")
    ag = AnnGeno(filename=anngeno_path, filemode="r", low_mem=True)
    
    print(f"Filtering for variants with MAF < {maf}")
    variants_to_keep_df = ag.annotations.filter((pl.col('AF_ukb') < maf))
    ag.subset_variants(set(variants_to_keep_df.select(pl.col("id")).collect()['id']))

    if only_snps:
        print(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) & 
            (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(snp_variants.select(pl.col('id')))

    associations_df = pl.read_parquet(associations_df_path)
    if debug:
        print("Debug is True, using only 5 associations")
        associations_df = associations_df.head()
    genes = associations_df['gene_id'].unique()

    gene_id_list = list(genes)
    valid_genes = [g for g in gene_id_list if g in ag.region_ids]
    invalid_regions = [g for g in gene_id_list if g not in ag.region_ids]
    if invalid_regions:
        print(f"Regions not found. Skipping regions {invalid_regions}")

    # Get regions in batches
    results = []
    for batch_genes in tqdm([valid_genes[i:i + batch_size] for i in range(0, len(valid_genes), batch_size)], desc="Loading region batches"):
        regions_dict = ag.get_many_regions(batch_genes)
        
        if device == "cuda":
            print("Using CUDA for computations.")
            batch_results = [get_gene_burdens_torch(regions_dict[gene]['genotypes'], regions_dict[gene]['annotations'], annotation_list, max_burden) for gene in tqdm(batch_genes, desc="Getting gene burdens")]

        else:
            print("Using CPU for computations.")
            batch_results = [get_gene_burdens_numba(regions_dict[gene]['genotypes'], regions_dict[gene]['annotations'], annotation_list, max_burden) for gene in tqdm(batch_genes, desc="Getting gene burdens")]
            # batch_results = Parallel(n_jobs=n_jobs, verbose=10)(
            #     delayed(get_gene_burdens)(regions_dict[gene]['genotypes'], regions_dict[gene]['annotations'], annotation_list, max_burden)
            #     for gene in tqdm(batch_genes) # tqdm for overall progress
            # )
        
        results.extend(batch_results)

    gene_burdens_sum = []
    gene_burdens_max = []
    gene_burdens_top2 = []

    for gis_sum, gis_max, gis_top2 in results:
        gene_burdens_sum.append(gis_sum)
        if max_burden:
            gene_burdens_max.append(gis_max)
            gene_burdens_top2.append(gis_top2)

    gene_burdens_sum_df = np.stack(gene_burdens_sum, axis=1)  # (n_samples, n_genes, n_annotations)
    if max_burden:
        gene_burdens_max_df = np.stack(gene_burdens_max, axis=1)  # (n_samples, n_genes, n_annotations)
        gene_burdens_top2_df = np.stack(gene_burdens_top2, axis=1)  # (n_samples, n_genes, n_annotations)
        print("Returning burdens: sum, max, and sum of top2")
        return gene_burdens_sum_df, gene_burdens_max_df, gene_burdens_top2_df, ag.samples, gene_id_list

    print("Max burden is False, returning only sum burden.")
    return gene_burdens_sum_df, gene_burdens_sum_df, gene_burdens_sum_df, ag.samples, gene_id_list


def compute_and_store_burdens(
    config_path,
    output_zarr,
    anno_scores_path=None,
    max_burden=False,
    only_snps=False,
    overwrite=False,
    debug=False,
    batch_size=32,
    n_jobs=32,
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
        max_burden (bool, optional): Whether to compute and store the maximum and sum(top2) burdens.
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

    if os.path.exists(zarr_file_path) and overwrite:
        print(f"Overwriting existing Zarr file at {zarr_file_path}.")
        try:
            shutil.rmtree(zarr_file_path)
        except Exception as e:
            print(f"Error deleting existing Zarr file: {e}")
            sys.exit(1)

    if os.path.exists(zarr_file_path):
        print(f"Zarr file exists at {zarr_file_path}, checking for new annotations.")
        try:
            root = zarr.group(zarr_file_path, mode='r+')
            sum_burdens = root.get("sum_burdens")
            zarr_annotations = root.get("annotations")
            max_burdens = root.get("max_burdens")
            top2_burdens = root.get("top2_burdens")

            if sum_burdens is None or zarr_annotations is None:
                print("Existing Zarr file is incomplete. Consider overwriting or creating a new one.")
                sys.exit(1)

            existing_annotations = list(zarr_annotations[:])
            new_annotation_list = list(set(all_annotation_list) - set(existing_annotations))
            union_annotation_list = existing_annotations + new_annotation_list

            if new_annotation_list:
                print(f"New annotations found: {new_annotation_list}")
                get_burdens_kwargs = {
                    "anngeno_path": config.get("anngeno_file"),
                    "associations_df_path": associations_df_path,
                    "maf": config.get("maf_upper_bound"),
                    "annotation_list": new_annotation_list,
                    "max_burden": max_burden,
                    "only_snps": only_snps,
                    "debug": debug,
                    "batch_size": batch_size,
                    "n_jobs": n_jobs,
                }
                if anno_scores_path:
                    get_burdens_kwargs["new_anno_df"] = anno_scores_df

                new_gene_burdens_sum_df, new_gene_burdens_max_df, new_gene_burdens_top2_df, _, _ = get_burdens_array(**get_burdens_kwargs)

                current_shape = sum_burdens.shape
                new_shape = (current_shape[0], current_shape[1], current_shape[2] + len(new_annotation_list))
                sum_burdens.resize(new_shape)
                sum_burdens[:, :, current_shape[2]:] = new_gene_burdens_sum_df

                if max_burden and max_burdens is not None:
                    print("Max burden is True, appending max_burdens array.")
                    max_burdens.resize(new_shape)
                    max_burdens[:, :, current_shape[2]:] = new_gene_burdens_max_df
                
                if max_burden and top2_burdens is not None:
                    print("Max burden is True, appending top2_burdens array.")
                    top2_burdens.resize(new_shape)
                    top2_burdens[:, :, current_shape[2]:] = new_gene_burdens_top2_df

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
            "batch_size": batch_size,
            "n_jobs": n_jobs,
        }
        if anno_scores_path:
            get_burdens_kwargs["new_anno_df"] = anno_scores_df

        gene_burdens_sum_df, gene_burdens_max_df, gene_burdens_top2_df, sample_id_arr, gene_id_list = get_burdens_array(**get_burdens_kwargs)

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

                root.create_array(
                    "top2_burdens",
                    shape=gene_burdens_top2_df.shape,
                    dtype=gene_burdens_top2_df.dtype,
                    chunks=(gene_burdens_top2_df.shape[0], gene_burdens_top2_df.shape[1], anno_chunk_size),
                    overwrite=overwrite,
                )[:] = gene_burdens_top2_df

            root.create_array("samples", shape=sample_id_arr.shape, dtype="str", overwrite=overwrite)[:] = sample_id_arr
            root.create_array("genes", shape=len(gene_id_list), dtype="str", overwrite=overwrite)[:] = gene_id_list
            root.create_array("annotations", shape=len(all_annotation_list), dtype="str", overwrite=overwrite)[:] = all_annotation_list
            
            print(f"Created new zarr array sum_burdens, {'max_burdens,' if max_burden else ''} samples, genes, and annotations.")

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
@click.option("--debug", is_flag=True, default=False, help="Use only 5 associations to debug code")
def compute_burdens(
    config_path: str,
    output_zarr: str,
    max_burden: bool = False,
    only_snps: bool = False,
    overwrite: bool = False,
    debug: bool = False,
    batch_size: int = 32,
    n_jobs: int = 32,
):
    print('You are running the script to compute gene burdens')

    compute_and_store_burdens(
        config_path=config_path,
        output_zarr=output_zarr,
        max_burden=max_burden,
        only_snps=only_snps,
        overwrite=overwrite,
        debug=debug,
        batch_size=batch_size,
        n_jobs=n_jobs,
    )

    print('Gene burdens have been computed and stored')

if __name__ == "__main__":
    cli()