import os
import gc
import sys
import yaml
import polars as pl
import numpy as np

import pathlib
from datetime import datetime

from tqdm import tqdm
from anngeno import AnnGeno

from numba import njit, prange
import multiprocessing

import logging

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)

logger = logging.getLogger(__name__)


@njit(parallel=True)
def compute_max_and_top2_chunked(
    score_vec, region_genotypes, chunk_size, no_variant_mask, valid_variant_idx
):
    """
    Compute per-sample max and sum of top 2 burdens for a single annotation,
    handling multiple variant copies correctly.
    """
    n_samples = region_genotypes.shape[1]
    n_chunks = (n_samples + chunk_size - 1) // chunk_size

    max_vals = np.empty(n_samples, dtype=np.float32)
    top2_sums = np.empty(n_samples, dtype=np.float32)

    for c in prange(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, n_samples)

        for s in range(start, end):
            if no_variant_mask[s]:
                max_vals[s] = np.nan
                top2_sums[s] = np.nan
                continue

            # Allocate buffer for this sample only
            max_expanded_size = (
                len(valid_variant_idx) * 2
            )  # at most 2 copies per variant
            expanded = np.zeros(max_expanded_size, dtype=np.float32)

            # Expand the burdens for this sample to account for homozygous variants
            idx_exp = 0
            for idx in valid_variant_idx:
                copies = int(region_genotypes[idx, s])
                if copies > 0:
                    burden = abs(score_vec[idx])
                    for _ in range(copies):
                        expanded[idx_exp] = burden
                        idx_exp += 1

            if idx_exp == 0:
                max_vals[s] = 0.0
                top2_sums[s] = 0.0
            elif idx_exp == 1:
                max_vals[s] = expanded[0]
                top2_sums[s] = expanded[0]
            else:
                top2 = np.partition(expanded, -2)[
                    -2:
                ]  # TODO: maybe we need np.partition(expanded[:idx_exp], -2)
                max_vals[s] = top2.max()
                top2_sums[s] = top2.sum()

    return max_vals, top2_sums


def get_gene_burdens_numba(
    region_genotypes,
    region_annotations,
    annotation_list,
    max_burden=False,
    chunk_size=None,
    na_mask=False,
):
    """
    Compute sum, max, and sum-of-top-2 burdens for all annotations in a gene region.
    Uses numba for speed.
    """
    # Identify samples with no variants if needed
    if na_mask:
        no_variant_mask = region_genotypes.sum(axis=0) == 0  # (samples, )
    else:
        no_variant_mask = np.zeros(region_genotypes.shape[1], dtype=bool)

    # Load and format variant scores (annotations)
    try:
        var_scores = (
            region_annotations[annotation_list]
            .fill_nan(0)
            .to_numpy()
            .astype(np.float32)
            .T
        )
    except Exception as e:
        print(f"Error: {e}\nReturning NaNs.")
        return np.nan, np.nan, np.nan

    # Sum burden: simple matrix multiply
    gis_sum = np.dot(var_scores, region_genotypes)
    if na_mask:
        gis_sum[:, no_variant_mask] = np.nan

    if not max_burden:
        return gis_sum, np.nan, np.nan

    # Chunking for parallel execution
    if chunk_size is None:
        num_cores = multiprocessing.cpu_count()
        chunk_size = max(1, region_genotypes.shape[1] // (num_cores * 2))

    # Find non-zero variants to skip unnecessary work
    nonzero_variants = np.any(region_genotypes > 0, axis=1)
    valid_variant_idx = np.where(nonzero_variants)[0].astype(np.int32)

    # Compute max + top2 for each annotation
    gis_max_list = []
    gis_top2_sum_list = []
    for a in tqdm(
        range(var_scores.shape[0]),
        desc=f"Computing max and top2, {var_scores.shape[1]} variants",
    ):
        score_vec = var_scores[a, :]
        max_vals, top2_sum = compute_max_and_top2_chunked(
            score_vec, region_genotypes, chunk_size, no_variant_mask, valid_variant_idx
        )
        gis_max_list.append(max_vals)
        gis_top2_sum_list.append(top2_sum)

        del max_vals, top2_sum
        gc.collect()

    # Stack all annotations into final arrays
    gis_max = np.stack(gis_max_list, axis=0)
    gis_top2 = np.stack(gis_top2_sum_list, axis=0)

    if na_mask:
        gis_max[:, no_variant_mask] = np.nan
        gis_top2[:, no_variant_mask] = np.nan

    return gis_sum, gis_max, gis_top2


def get_burdens_array_streaming(
    anngeno, gene_id_list, annotation_list, max_burden=False, device="cuda"
):
    with torch.no_grad():
        G = torch.tensor(
            region_genotypes, dtype=torch.float32, device=device
        )  # (samples, variants)
        A = torch.tensor(
            region_annotations[annotation_list].fill_nan(0).to_numpy(),
            dtype=torch.float32,
            device=device,
        )  # (variants, annotations)

        no_variant_mask = G.sum(dim=1) == 0
        gis_sum = G @ A
        gis_sum[no_variant_mask] = float("nan")

        if not max_burden:
            return gis_sum.cpu().numpy(), None, None

        gis_max = []
        gis_top2 = []

        for a in range(A.shape[1]):
            scores = torch.abs(G * A[:, a])  # (samples, variants)
            scores[G == 0] = float("-inf")

            top2_vals, _ = torch.topk(scores, k=2, dim=1, largest=True, sorted=False)
            gis_max.append(torch.max(top2_vals, dim=1).values)
            gis_top2.append(top2_vals.sum(dim=1))

            torch.cuda.empty_cache()

        gis_max = torch.stack(gis_max, dim=1)
        gis_top2 = torch.stack(gis_top2, dim=1)

        gis_max[no_variant_mask] = float("nan")
        gis_top2[no_variant_mask] = float("nan")

        return gis_sum.cpu().numpy(), gis_max.cpu().numpy(), gis_top2.cpu().numpy()


def get_burdens_array(
    config,
    associations_df_path,
    annotation_list,
    new_anno_df=None,
    max_burden=False,
    only_snps=False,
    debug=False,
    n_jobs=-1,
    batch_size=100,
    device="cuda" if torch.cuda.is_available() else "cpu",
):
    maf = config.get("maf_upper_bound")

    print("Loading AnnGeno file")
    anngeno_file = config.get("anngeno_file")
    ag = AnnGeno(filename=anngeno_file, filemode="r", low_mem=True)

    print(f"Filtering for variants with MAF < {maf}")

    variants_to_keep_df = ag.annotations.filter((pl.col("MAF") < maf))
    ag.subset_variants(set(variants_to_keep_df.select(pl.col("id")).collect()["id"]))

    if only_snps:
        print(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(snp_variants.select(pl.col("id")))

    associations_df = pl.read_parquet(associations_df_path)
    if debug:
        print("Debug is True, using only 5 associations")
        associations_df = associations_df.head()
    genes = associations_df["region"].unique()

    gene_id_list = list(genes)
    valid_genes = [g for g in gene_id_list if g in ag.region_ids]
    invalid_regions = [g for g in gene_id_list if g not in ag.region_ids]
    if invalid_regions:
        print(f"Regions not found. Skipping regions {invalid_regions}")

    # Get regions in batches
    results = []
    for batch_genes in tqdm(
        [
            valid_genes[i : i + batch_size]
            for i in range(0, len(valid_genes), batch_size)
        ],
        desc="Loading region batches",
    ):
        regions_dict = ag.get_many_regions(batch_genes)

        if device == "cuda":
            print("Using CUDA for computations.")
            batch_results = [
                get_gene_burdens_torch(
                    regions_dict[gene]["genotypes"],
                    regions_dict[gene]["annotations"],
                    annotation_list,
                    max_burden,
                )
                for gene in tqdm(batch_genes, desc="Getting gene burdens")
            ]

        else:
            print("Using CPU for computations.")
            results = Parallel(n_jobs=n_jobs, verbose=10)(
                delayed(get_gene_burdens)(
                    regions_dict[gene]["genotypes"],
                    regions_dict[gene]["annotations"],
                    annotation_list,
                    max_burden,
                )
                for gene in tqdm(batch_genes)  # tqdm for overall progress
            )

        results.extend(batch_results)

    gene_burdens_sum = []
    gene_burdens_max = []
    gene_burdens_top2 = []

    for gis_sum, gis_max, gis_top2 in results:
        gene_burdens_sum.append(gis_sum)
        if max_burden:
            gene_burdens_max.append(gis_max)
            gene_burdens_top2.append(gis_top2)

    gene_burdens_sum_df = np.stack(
        gene_burdens_sum, axis=1
    )  # (n_samples, n_genes, n_annotations)
    if max_burden:
        gene_burdens_max_df = np.stack(
            gene_burdens_max, axis=1
        )  # (n_samples, n_genes, n_annotations)
        gene_burdens_top2_df = np.stack(
            gene_burdens_top2, axis=1
        )  # (n_samples, n_genes, n_annotations)
        print("Returning burdens: sum, max, and sum of top2")
        return (
            gene_burdens_sum_df,
            gene_burdens_max_df,
            gene_burdens_top2_df,
            ag.samples,
            gene_id_list,
        )

    print("Max burden is False, returning only sum burden.")
    return (
        gene_burdens_sum_df,
        gene_burdens_sum_df,
        gene_burdens_sum_df,
        ag.samples,
        gene_id_list,
    )


def compute_and_store_burdens(
    config_path,
    gene_list,
    output_dir,
    new_annotation_df=None,
    only_snps=False,
    variant_subset=None,
    sample_subset=None,
    max_burden=True,
    gene_chunk_size=50,
    sample_chunk_size=5_000,
    na_mask=False,
    overwrite=False,
    debug=False,
):
    """
    Computes and stores gene burdens directly to Zarr in streaming mode.
    Either creates new Zarr arrays or overwrites existing ones.
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    logger.info("Loading AnnGeno file")
    ag = AnnGeno(filename=config.get("anngeno_file"), filemode="r", low_mem=True)

    if new_annotation_df is not None:
        logger.info("Setting new annotations to AnnGeno")
        if isinstance(new_annotation_df, pl.LazyFrame):
            ag._set_annotations(new_annotation_df)
        else:
            ag._set_annotations(new_annotation_df.lazy())

    if variant_subset:
        logger.info(
            f"Filtering for variants in subset of {len(variant_subset)} variants"
        )
        variants_to_keep = (
            ag.annotations.filter(pl.col("id").is_in(set(variant_subset)))
            .select("id")
            .collect()["id"]
        )
        ag.subset_variants(set(variants_to_keep))

    maf = config.get("maf", None)
    if maf is not None:
        logger.info(f"Filtering for variants with MAF < {maf}")
        variants_to_keep = (
            ag.annotations.filter((pl.col("AF_ukb") < maf)).select("id").collect()["id"]
        )
        ag.subset_variants(set(variants_to_keep))

    if only_snps:
        logger.info(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(set(snp_variants.select(pl.col("id")).collect()["id"]))

    if sample_subset:
        logger.info(
            f"Filtering for samples. Restricting to {len(sample_subset)} samples"
        )
        ag.subset_samples(sample_subset)

    logger.info("Drop is_nans from annotations")
    sel_cols = [
        col
        for col in ag.annotations.collect_schema().names()
        if not col.endswith("is_nan")
    ]
    ag.subset_annotations(sel_cols)

    all_annotation_list = []

    if anno_scores_path:
        anno_scores_df = pl.read_parquet(anno_scores_path)
        available_annotations = list(
            set(anno_scores_df.columns) - set(["chrom", "pos", "ref", "alt", "region"])
        )
        all_annotation_list.extend(available_annotations)
    else:
        rare_variant_annotations_dict = config.get("rare_variant_annotations")
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
            root = zarr.group(zarr_file_path, mode="r+")
            sum_burdens = root.get("sum_burdens")
            zarr_annotations = root.get("annotations")
            max_burdens = root.get("max_burdens")
            top2_burdens = root.get("top2_burdens")

            if sum_burdens is None or zarr_annotations is None:
                print(
                    "Existing Zarr file is incomplete. Consider overwriting or creating a new one."
                )
                sys.exit(1)

            existing_annotations = list(zarr_annotations[:])
            new_annotation_list = list(
                set(all_annotation_list) - set(existing_annotations)
            )
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

                (
                    new_gene_burdens_sum_df,
                    new_gene_burdens_max_df,
                    new_gene_burdens_top2_df,
                    _,
                    _,
                ) = get_burdens_array(**get_burdens_kwargs)

                current_shape = sum_burdens.shape
                new_shape = (
                    current_shape[0],
                    current_shape[1],
                    current_shape[2] + len(new_annotation_list),
                )
                sum_burdens.resize(new_shape)
                sum_burdens[:, :, current_shape[2] :] = new_gene_burdens_sum_df

                if max_burden and max_burdens is not None:
                    print("Max burden is True, appending max_burdens array.")
                    max_burdens.resize(new_shape)
                    max_burdens[:, :, current_shape[2] :] = new_gene_burdens_max_df

                if max_burden and top2_burdens is not None:
                    print("Max burden is True, appending top2_burdens array.")
                    top2_burdens.resize(new_shape)
                    top2_burdens[:, :, current_shape[2] :] = new_gene_burdens_top2_df

                if zarr_annotations is not None:
                    zarr_annotations_new_shape = (len(union_annotation_list),)
                    if zarr_annotations.shape != zarr_annotations_new_shape:
                        zarr_annotations.resize(zarr_annotations_new_shape)
                    zarr_annotations[:] = union_annotation_list
                    print("Updated 'annotations' array.")
                print(f"Appended data for new annotations: {new_annotation_list}")
            else:
                # File doesn't exist; safe to write
                df_lazy.sink_parquet(gene_file)

        except Exception as e:
            print(f"Error computing burdens: {e}")
            sys.exit(1)

        gc.collect()

    logger.info(
        f"Stored burdens for {n_genes} genes and {n_annos} annotations in {output_dir}"
    )


import click


@click.command()
@click.option(
    "--config-path",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML config file.",
)
@click.option(
    "--gene-list-path",
    required=True,
    type=click.Path(exists=True),
    help="List of genes to compute the burdens for.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Directory to write per-gene Parquet files.",
)
@click.option("--only-snps", is_flag=True, default=False, help="Filter for SNPs only.")
@click.option(
    "--variant-subset-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to text file with variant IDs to include.",
)
@click.option(
    "--sample-subset-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to text file with sample IDs to include.",
)
@click.option(
    "--gene-chunk-size",
    type=int,
    default=50,
    help="Number of genes to process per chunk.",
)
@click.option(
    "--sample-chunk-size",
    type=int,
    default=5000,
    help="Number of samples to process per chunk.",
)
@click.option(
    "--na-mask",
    is_flag=True,
    default=False,
    help="Filter out samples with no variants in the region.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Whether to overwrite existing gene Parquet files.",
)
def cli(
    config_path,
    gene_list_path,
    output_dir,
    only_snps,
    variant_subset_path,
    sample_subset_path,
    gene_chunk_size,
    sample_chunk_size,
    na_mask,
    overwrite,
):
    try:
        ext = pathlib.Path(gene_list_path).suffix.lower()
        if (ext == ".parquet") or (ext == ".pq"):
            gene_list = pl.read_parquet(gene_list_path)["gene_id"].unique().to_list()
        else:
            gene_list = (
                pl.read_csv(gene_list_path, has_header=False)
                .to_series()
                .unique()
                .to_list()
            )
    except Exception as e:
        logger.error(f"No gene list provided or error reading gene list: {e}")
        sys.exit(1)

    if variant_subset_path:
        variant_subset = (
            pl.read_csv(variant_subset_path, has_header=False).to_series().to_list()
        )
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

        (
            gene_burdens_sum_df,
            gene_burdens_max_df,
            gene_burdens_top2_df,
            sample_id_arr,
            gene_id_list,
        ) = get_burdens_array(**get_burdens_kwargs)

        try:
            root = zarr.group(zarr_file_path)
            root.create_array(
                "sum_burdens",
                shape=gene_burdens_sum_df.shape,
                dtype=gene_burdens_sum_df.dtype,
                chunks=(
                    gene_burdens_sum_df.shape[0],
                    gene_burdens_sum_df.shape[1],
                    anno_chunk_size,
                ),
                overwrite=overwrite,
            )[:] = gene_burdens_sum_df

            if max_burden:
                root.create_array(
                    "max_burdens",
                    shape=gene_burdens_max_df.shape,
                    dtype=gene_burdens_max_df.dtype,
                    chunks=(
                        gene_burdens_max_df.shape[0],
                        gene_burdens_max_df.shape[1],
                        anno_chunk_size,
                    ),
                    overwrite=overwrite,
                )[:] = gene_burdens_max_df

                root.create_array(
                    "top2_burdens",
                    shape=gene_burdens_top2_df.shape,
                    dtype=gene_burdens_top2_df.dtype,
                    chunks=(
                        gene_burdens_top2_df.shape[0],
                        gene_burdens_top2_df.shape[1],
                        anno_chunk_size,
                    ),
                    overwrite=overwrite,
                )[:] = gene_burdens_top2_df

            root.create_array(
                "samples", shape=sample_id_arr.shape, dtype="str", overwrite=overwrite
            )[:] = sample_id_arr
            root.create_array(
                "genes", shape=len(gene_id_list), dtype="str", overwrite=overwrite
            )[:] = gene_id_list
            root.create_array(
                "annotations",
                shape=len(all_annotation_list),
                dtype="str",
                overwrite=overwrite,
            )[:] = all_annotation_list

            print(
                f"Created new zarr array sum_burdens, {'max_burdens,' if max_burden else ''} samples, genes, and annotations."
            )

        except Exception as e:
            print(f"Error creating new Zarr file: {e}")
            sys.exit(1)


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--config-path", type=str, required=True, help="Config file with all details"
)
@click.option("--output-zarr", type=str, required=True, help="Output zarr file")
@click.option(
    "--max-burden", is_flag=True, default=False, help="Compute burdens using max scores"
)
@click.option(
    "--only-snps", is_flag=True, default=False, help="Use only SNPs to compute burdens"
)
@click.option(
    "--overwrite", is_flag=True, default=False, help="Overwrite old zarr file"
)  # TODO
@click.option(
    "--debug", is_flag=True, default=False, help="Use only 5 associations to debug code"
)
def compute_burdens(
    config_path: str,
    output_zarr: str,
    max_burden: bool = False,
    only_snps: bool = False,
    overwrite: bool = False,
    debug: bool = False,
):
    print("You are running the script to compute gene burdens")

    compute_and_store_burdens(
        config_path=config_path,
        gene_list=gene_list,
        output_dir=output_dir,
        only_snps=only_snps,
        overwrite=overwrite,
        debug=debug,
    )


if __name__ == "__main__":
    cli()
