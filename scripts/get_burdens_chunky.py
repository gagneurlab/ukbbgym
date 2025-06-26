import gc
import sys
import yaml
import zarr
import click
import polars as pl
import numpy as np

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
def compute_max_and_top2_chunked(score_vec, region_genotypes, chunk_size, no_variant_mask):
    n_variants, n_samples = region_genotypes.shape
    n_chunks = (n_samples + chunk_size - 1) // chunk_size

    max_vals = np.empty(n_samples, dtype=np.float32)
    top2_sums = np.empty(n_samples, dtype=np.float32)

    for c in prange(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, n_samples)
        for s in range(start, end):
            if no_variant_mask[s]:
                    # Skip computation, set NaN
                    max_vals[s] = np.nan
                    top2_sums[s] = np.nan
                    continue
            
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
    no_variant_mask = region_genotypes.sum(axis = 0) == 0 # (samples, )

    try:
        var_scores = region_annotations[annotation_list].fill_nan(0).to_numpy().astype(np.float32).transpose()  # shape: (annotations, variants)
    except Exception as e:
        logger.info(f"Error: {e}\nReturning NaNs.")
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

    # Compute max + top2 via numba
    gis_max_list = []
    gis_top2_sum_list = []
    for a in range(var_scores.shape[0]):
        score_vec = var_scores[a, :]
        max_vals, top2_sum = compute_max_and_top2_chunked(score_vec, region_genotypes, chunk_size, no_variant_mask)
        gis_max_list.append(max_vals)
        gis_top2_sum_list.append(top2_sum)
        
        del max_vals, top2_sum
        gc.collect()

    gis_max = np.stack(gis_max_list, axis=0).T         # (samples, annotations)
    gis_top2 = np.stack(gis_top2_sum_list, axis=0).T   # (samples, annotations)

    gis_max[no_variant_mask, :] = np.nan
    gis_top2[no_variant_mask, :] = np.nan

    return gis_sum, gis_max, gis_top2


def get_burdens_array_streaming(
    anngeno,
    gene_id_list,
    annotation_list,
    gene_chunk_size=50,
    sample_slice=None,
    max_burden=True,  #TODO
    device="cpu",     #TODO
):
    """
    Generator yielding (gene, sum_burden, max_burden, top2_burden) for each gene.
    """

    for i in tqdm(range(0, len(gene_id_list), gene_chunk_size), position=0, leave=True, desc=f"Genes: Processing chunks of {gene_chunk_size} genes"):
        batch_genes = gene_id_list[i : i + gene_chunk_size]
        regions_dict = anngeno.get_many_regions(
            regions=batch_genes, 
            sample_slice=sample_slice,
            )

        for gene in batch_genes:
            burdens = get_gene_burdens_numba(
                regions_dict[gene]["genotypes"],
                regions_dict[gene]["annotations"],
                annotation_list=annotation_list,
                max_burden=True,
            )
            yield gene, *burdens
            del burdens
            gc.collect()
        del regions_dict
        gc.collect()

def compute_and_store_burdens(
    config_path,
    associations_df_path,
    output_zarr,
    only_snps=False,
    sample_set=None,
    gene_chunk_size=50,
    sample_chunk_size=5_000,
    device="cpu",
    overwrite=True,
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

    if only_snps:
        logger.info(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) &
            (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(snp_variants.select(pl.col('id')))

    if sample_set:
        logger.info(f"Filtering for samples. Restricting to {len(sample_set)} samples")
        ag.subset_samples(sample_set)

    all_annotation_list = []
    rare_variant_annotations_dict = config.get('rare_variant_annotations')
    if rare_variant_annotations_dict:
        for category in rare_variant_annotations_dict.values():
            all_annotation_list.extend(category)
    all_annotation_list = list(set(all_annotation_list).intersection(set(ag.annotations.collect_schema().names())))

    associations_df = pl.read_parquet(associations_df_path)
    gene_id_list = associations_df['gene_id'].unique()
    valid_genes = [g for g in gene_id_list if g in ag.region_ids]
    invalid_regions = [g for g in gene_id_list if g not in ag.region_ids]
    if invalid_regions:
        logger.info(f"Regions not found. Skipping regions {invalid_regions}")

    sample_ids = ag.samples
    n_samples = len(sample_ids)
    n_genes = len(valid_genes)
    n_annos = len(all_annotation_list)

    mode = "w" if overwrite else "a"
    zarr_root = zarr.open_group(output_zarr, mode=mode)

    logger.info(f"Creating new Zarr arrays at {output_zarr}")
    sum_burdens = zarr_root.create_array("sum_burdens", shape=(n_samples, n_genes, n_annos), chunks=(sample_chunk_size, 1, 1), dtype="f4")
    max_burdens = zarr_root.create_array("max_burdens", shape=(n_samples, n_genes, n_annos), chunks=(sample_chunk_size, 1, 1), dtype="f4")
    top2_burdens = zarr_root.create_array("top2_burdens", shape=(n_samples, n_genes, n_annos), chunks=(sample_chunk_size, 1, 1), dtype="f4")
    zarr_root.create_array("samples", shape=(n_samples,), dtype="U50")
    zarr_root.create_array("annotations", shape=(n_annos,), chunks=(n_annos,), dtype="U50")
    zarr_root.create_array("genes", shape=(n_genes,), chunks=(1,), dtype="U50")

    zarr_root["annotations"][:] = np.array(all_annotation_list, dtype="U50")
    zarr_root["samples"][:] = sample_ids
    gene_array = zarr_root["genes"]
    gene_idx_map = {g: i for i, g in enumerate(valid_genes)}

    for start in range(0, n_samples, sample_chunk_size):
        end = min(start + sample_chunk_size, n_samples)
        sample_slice = slice(start, end)

        logger.info(f"Processing samples {start}:{end} ({end-start} samples)")

        for gene, s_burden, m_burden, t2_burden in get_burdens_array_streaming(
            ag,
            valid_genes,
            all_annotation_list,
            gene_chunk_size=gene_chunk_size,
            sample_slice=sample_slice,
            device=device,
        ):
            idx_gene = gene_idx_map[gene]
            sum_burdens[sample_slice, idx_gene, :] = s_burden
            max_burdens[sample_slice, idx_gene, :] = m_burden
            top2_burdens[sample_slice, idx_gene, :] = t2_burden
            gene_array[idx_gene] = gene

        gc.collect()

    logger.info(f"Stored burdens for {n_genes} genes and {n_annos} annotations in Zarr at {output_zarr}")
