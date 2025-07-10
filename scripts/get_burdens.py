import os
import gc
import sys
import yaml
import polars as pl
import numpy as np

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
def compute_max_and_top2_chunked(score_vec, region_genotypes, chunk_size, no_variant_mask, valid_variant_idx):
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
            max_expanded_size = len(valid_variant_idx) * 2  # at most 2 copies per variant
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
                top2 = np.partition(expanded, -2)[-2:] #TODO: maybe we need np.partition(expanded[:idx_exp], -2)
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
        var_scores = region_annotations[annotation_list].fill_nan(0).to_numpy().astype(np.float32).T
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
    for a in tqdm(range(var_scores.shape[0]), desc=f"Computing max and top2, {var_scores.shape[1]} variants"):
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
    anngeno,
    gene_id_list,
    annotation_list,
    gene_chunk_size=50,
    sample_slice=None,
    na_mask=False,
    max_burden=True,  #TODO
    device="cpu",     #TODO
):
    """
    Generator yielding (gene, sum_burden, max_burden, top2_burden) for each gene.
    """

    for i in range(0, len(gene_id_list), gene_chunk_size):
        batch_genes = gene_id_list[i : i + gene_chunk_size]
        print(f"{[datetime.now().strftime('%Y-%m-%d %H:%M:%S')]} Starting loading {gene_chunk_size} regions")
        regions_dict = anngeno.get_many_regions(
            regions=batch_genes, 
            sample_slice=sample_slice,
            observed_only=True,
            )
        print(f"{[datetime.now().strftime('%Y-%m-%d %H:%M:%S')]} Done loading {gene_chunk_size} regions")

        for gene in batch_genes:
            burdens = get_gene_burdens_numba(
                regions_dict[gene]["genotypes"],
                regions_dict[gene]["annotations"],
                annotation_list=annotation_list,
                max_burden=True,
                na_mask=na_mask,
            )
            yield gene, *burdens
            del burdens
            gc.collect()
        del regions_dict
        gc.collect()


def compute_and_store_burdens(
    config_path,
    gene_list,
    output_dir,
    only_snps=False,
    sample_set=None,
    gene_chunk_size=50,
    sample_chunk_size=5_000,
    na_mask=False,
    device="cpu",
    overwrite=False,
):
    """
    Computes and stores gene burdens directly to Zarr in streaming mode.
    Either creates new Zarr arrays or overwrites existing ones.
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    logger.info("Loading AnnGeno file")
    ag = AnnGeno(filename=config.get("anngeno_file"), filemode="r", low_mem=True)
    maf = config.get('maf', 0.001)

    logger.info(f"Filtering for variants with MAF < {maf}")
    variants_to_keep_df = ag.annotations.filter((pl.col('AF_ukb') < maf))
    ag.subset_variants(set(variants_to_keep_df.select(pl.col("id")).collect()['id']))
    
    logger.info("Drop is_nans from annotations")
    sel_cols = [col for col in ag.annotations.collect_schema().names() if not col.endswith('is_nan')]
    ag.subset_annotations(sel_cols)

    if only_snps:
        logger.info(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) &
            (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(set(snp_variants.select(pl.col('id')).collect()['id']))

    if sample_set:
        logger.info(f"Filtering for samples. Restricting to {len(sample_set)} samples")
        ag.subset_samples(sample_set)

    all_annotation_list = []
    rare_variant_annotations_dict = config.get('rare_variant_annotations')
    if rare_variant_annotations_dict:
        for category in rare_variant_annotations_dict.values():
            all_annotation_list.extend(category)
    all_annotation_list = list(set(all_annotation_list).intersection(set(ag.annotations.collect_schema().names())))

    # Get valid genes from the provided gene list
    valid_genes = [g for g in gene_list if g in ag.region_ids]
    invalid_regions = [g for g in gene_list if g not in ag.region_ids]
    if invalid_regions:
        logger.info(f"Regions not found. Skipping regions {invalid_regions}")

    sample_ids = ag.samples
    n_samples = len(sample_ids)
    n_genes = len(valid_genes)
    n_annos = len(all_annotation_list)

    logger.info(f"Found {n_genes} valid genes and {n_annos} annotations across {n_samples} samples.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files_overwritten = set() # To track files that have been overwritten

    for start in tqdm(range(0, n_samples, sample_chunk_size), desc=f"Processing {sample_chunk_size} sample chunks"):
        end = min(start + sample_chunk_size, n_samples)
        sample_slice = slice(start, end)
        current_sample_ids = ag.samples[sample_slice.start:sample_slice.stop]

        for gene, s_burden, m_burden, t2_burden in get_burdens_array_streaming(
            ag,
            valid_genes,
            all_annotation_list,
            gene_chunk_size=gene_chunk_size,
            sample_slice=sample_slice,
            na_mask=na_mask,
            device=device,
        ):
            n_samples_in_chunk = s_burden.shape[1]
            annotation_ids = np.array(all_annotation_list)

            sample_col = np.tile(current_sample_ids, n_annos)
            annotation_col = np.repeat(annotation_ids, n_samples_in_chunk)

            df_lazy = pl.LazyFrame({
                "sample_id": sample_col,
                "gene_id": [gene] * (n_samples_in_chunk * n_annos),
                "annotation": annotation_col,
                "sum": s_burden.flatten(),
                "max": m_burden.flatten(),
                "top2": t2_burden.flatten(),
            })

            gene_file = os.path.join(output_dir, f"{gene}.parquet")
            if os.path.exists(gene_file):
                if overwrite:
                    if gene_file not in files_overwritten:
                        # Overwrite only once if it existed before
                        df_lazy.sink_parquet(gene_file)
                        files_overwritten.add(gene_file)
                    else:
                        # File was already overwritten, so we append
                        try:
                            existing = pl.read_parquet(gene_file).lazy()
                            df_lazy = pl.concat([existing, df_lazy])
                            df_lazy.sink_parquet(gene_file)
                        except Exception as e:
                            logger.warning(f"Could not concat {gene}: {e}")
                else:
                    # No overwrite allowed, always append
                    try:
                        existing = pl.read_parquet(gene_file).lazy()
                        df_lazy = pl.concat([existing, df_lazy])
                        df_lazy.sink_parquet(gene_file)
                    except Exception as e:
                        logger.warning(f"Could not concat {gene}: {e}")
            else:
                # File doesn't exist; safe to write
                df_lazy.sink_parquet(gene_file)

        gc.collect()

    logger.info(f"Stored burdens for {n_genes} genes and {n_annos} annotations in {output_dir}")


import click
@click.command()
@click.option('--config-path', required=True, type=click.Path(exists=True), help="Path to YAML config file.")
@click.option('--gene-list', required=True, type=click.Path(exists=True), help="List of genes to compute the burdens for.")
@click.option('--output-dir', required=True, type=click.Path(), help="Directory to write per-gene Parquet files.")
@click.option('--only-snps', is_flag=True, default=False, help="Filter for SNPs only.")
@click.option('--sample-set-path', type=click.Path(exists=True), default=None, help="Optional path to text file with sample IDs to include.")
@click.option('--gene-chunk-size', type=int, default=50, help="Number of genes to process per chunk.")
@click.option('--sample-chunk-size', type=int, default=5000, help="Number of samples to process per chunk.")
@click.option('--na-mask', is_flag=True, default=False, help="Filter out samples with no variants in the region.")
@click.option('--device', default='cpu', help="Device to use for computation.")
@click.option('--overwrite', is_flag=True, default=False, help="Whether to overwrite existing gene Parquet files.")
def cli(
    config_path,
    gene_list,
    output_dir,
    only_snps,
    sample_set_path,
    gene_chunk_size,
    sample_chunk_size,
    na_mask,
    device,
    overwrite,
):
    if sample_set_path:
        with open(sample_set_path) as f:
            sample_set = [line.strip() for line in f if line.strip()]
    else:
        sample_set = None

    compute_and_store_burdens(
        config_path=config_path,
        gene_list=gene_list,
        output_dir=output_dir,
        only_snps=only_snps,
        sample_set=sample_set,
        gene_chunk_size=gene_chunk_size,
        sample_chunk_size=sample_chunk_size,
        na_mask=na_mask,
        device=device,
        overwrite=overwrite
    )

if __name__ == "__main__":
    cli()
