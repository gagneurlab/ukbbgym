import gc
import yaml
import zarr
import click
import polars as pl
import numpy as np

from tqdm import tqdm
from anngeno import AnnGeno

from numba import njit, prange
import multiprocessing

@njit(parallel=True)
def compute_max_and_top2_batch(score_matrix, genotype_tensor, chunk_size, no_variant_mask):
    n_genes, n_variants, n_samples = genotype_tensor.shape
    max_vals = np.empty((n_samples, n_genes), dtype=np.float32)
    top2_sums = np.empty((n_samples, n_genes), dtype=np.float32)

    for g in prange(n_genes):
        for s in range(n_samples):
            if no_variant_mask[g, s]:
                max_vals[s, g] = np.nan
                top2_sums[s, g] = np.nan
                continue
            burden = np.abs(score_matrix[g, :] * genotype_tensor[g, :, s])
            if burden.size >= 2:
                top2 = np.partition(burden, -2)[-2:]
                max_vals[s, g] = top2.max()
                top2_sums[s, g] = top2.sum()
            elif burden.size == 1:
                max_vals[s, g] = burden[0]
                top2_sums[s, g] = burden[0]
            else:
                max_vals[s, g] = 0.0
                top2_sums[s, g] = 0.0

    return max_vals, top2_sums

def get_gene_burdens_numba(
    region_genotypes_list,
    region_annotations_list,
    annotation_name,
):
    n_genes = len(region_genotypes_list)
    n_samples = region_genotypes_list[0].shape[1]

    genotype_tensor = np.stack(region_genotypes_list, axis=0)  # (genes, variants, samples)
    score_matrix = np.stack([
        region_annotations[annotation_name].fill_nan(0).to_numpy().astype(np.float32)
        for region_annotations in region_annotations_list
    ], axis=0)  # (genes, variants)

    no_variant_mask = np.stack([g.sum(axis=0) == 0 for g in region_genotypes_list], axis=0).T  # (samples, genes)

    # Sum burden
    sum_burdens = np.tensordot(score_matrix, genotype_tensor, axes=([1], [1])).transpose(2, 0).astype(np.float32)  # (samples, genes)
    sum_burdens[no_variant_mask] = np.nan

    max_burdens, top2_burdens = compute_max_and_top2_batch(score_matrix, genotype_tensor, 512, no_variant_mask)
    max_burdens[no_variant_mask] = np.nan
    top2_burdens[no_variant_mask] = np.nan

    return sum_burdens, max_burdens, top2_burdens  # (samples, genes)


def get_burdens_array_streaming(
    anngeno,
    gene_id_list,
    annotation_list,
    gene_chunk_size=50,
    sample_chunk_size=None,
    max_burden=True,  #TODO
    device="cpu",     #TODO
):
    """
    Generator yielding (gene, sum_burden, max_burden, top2_burden) for each gene.
    """
    
    # for i in tqdm(range(0, len(gene_id_list), gene_chunk_size)):
    #     batch_genes = gene_id_list[i : i + gene_chunk_size]
    #     regions_dict = anngeno.get_many_regions(
    #         regions=batch_genes, 
    #         sample_slice=sample_slice,
    #         )

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
    # del regions_dict
    # gc.collect()

def compute_and_store_burdens(
    config_path,
    associations_df_path,
    output_zarr,
    only_snps=False,
    sample_set=None,
    gene_chunk_size=50,
    sample_chunk_size=5_000,
    device="cpu",
):
    """
    Computes and stores gene burdens directly to Zarr in streaming mode,
    with extendable gene axis.
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print("Loading AnnGeno file")
    ag = AnnGeno(filename=config.get("anngeno_file"), filemode="r", low_mem=True)
    maf = config.get('maf', 0.001)
    
    print(f"Filtering for variants with MAF < {maf}")
    variants_to_keep_df = ag.annotations.filter((pl.col('AF_ukb') < maf))
    ag.subset_variants(set(variants_to_keep_df.select(pl.col("id")).collect ()['id']))
    
    if only_snps:
        print(f"Filtering for SNPs only")
        snp_variants = ag.annotations.filter(
            (pl.col("ref").str.len_chars() == 1) & 
            (pl.col("alt").str.len_chars() == 1)
        )
        ag.subset_variants(snp_variants.select(pl.col('id')))
    
    if sample_set:
        print(f"Filtering for samples. Restricting to {len(sample_set)} samples")
        ag.subset_samples(sample_set)

    all_annotation_list = []
    rare_variant_annotations_dict = config.get('rare_variant_annotations')
    if rare_variant_annotations_dict:
        for category in rare_variant_annotations_dict.values():
            all_annotation_list.extend(category)
    all_annotation_list = list(set(all_annotation_list).intersection(set(ag.annotations.collect_schema().names())))  # Ensure unique annotations are present in the AnnGeno object
    n_annos = len(all_annotation_list)

    associations_df = pl.read_parquet(associations_df_path)
    gene_id_list = associations_df['gene_id'].unique()
    valid_genes = [g for g in gene_id_list if g in ag.region_ids]
    invalid_regions = [g for g in gene_id_list if g not in ag.region_ids]
    if invalid_regions:
        print(f"Regions not found. Skipping regions {invalid_regions}")
    sample_ids = ag.samples
    n_samples = len(sample_ids)
    n_genes = len(valid_genes)

    if sample_chunk_size is None:
        sample_chunk_size = n_samples  # no chunking

    # compressors = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    zarr_root = zarr.open_group(output_zarr, mode="a")

    # ===========================
    # Initialize or Extend Zarr
    # ===========================

    if "sum_burdens" not in zarr_root:
        print(f"Creating new Zarr arrays at {output_zarr}")
        sum_burdens = zarr_root.create_array(
            "sum_burdens",
            shape=(n_samples, n_genes, 0),  # 0 annotations initially
            chunks=(sample_chunk_size, 1, 1),
            dtype="f4"
        )
        max_burdens = zarr_root.create_array(
            "max_burdens",
            shape=(n_samples, n_genes, 0),
            chunks=(sample_chunk_size, 1, 1),
            dtype="f4"
        )
        top2_burdens = zarr_root.create_array(
            "top2_burdens",
            shape=(n_samples, n_genes, 0),
            chunks=(sample_chunk_size, 1, 1),
            dtype="f4"
        )
        zarr_root.create_array("samples", shape=(n_samples,), dtype="U50")
        zarr_root.create_array("genes", shape=(n_genes,), dtype="U50")
        zarr_root.create_array("annotations", shape=(0,), dtype="U50")
        
        zarr_root["samples"][:] = np.array(sample_ids, dtype="U50")
        zarr_root["genes"][:] = np.array(valid_genes, dtype="U50")
    else:
        print(f"Reopening existing Zarr store: {output_zarr}")
        sum_burdens = zarr_root["sum_burdens"]
        max_burdens = zarr_root["max_burdens"]
        top2_burdens = zarr_root["top2_burdens"]

    # Handle annotation extension
    existing_annotations = list(zarr_root["annotations"][:])
    new_annotations = list(set(all_annotation_list) - set(existing_annotations))
    all_annotations_combined = existing_annotations + sorted(new_annotations)
    n_annos_total = len(all_annotations_combined)

    if len(new_annotations) > 0:
        print(f"Extending with {len(new_annotations)} new annotations: {new_annotations}")
        sum_burdens.resize((n_samples, n_genes, n_annos_total))
        max_burdens.resize((n_samples, n_genes, n_annos_total))
        top2_burdens.resize((n_samples, n_genes, n_annos_total))
        zarr_root["annotations"].resize((n_annos_total,))
        zarr_root["annotations"][:] = np.array(all_annotations_combined, dtype="U50")

    
    zarr_root["annotations"].resize(annotation_offset + n_annos)
    annotation_array = zarr_root["annotations"]
    annotation_idx_map = {a: i for i, a in enumerate(all_annotations_combined)}

    for anno in tqdm(new_annotations, desc="Annotations"):
    anno_idx = annotation_idx_map[anno]

    for g_start in range(0, len(valid_genes), gene_chunk_size):
        batch_genes = valid_genes[g_start : g_start + gene_chunk_size]
        gene_slice = slice(g_start, g_start + len(batch_genes))

        for s_start in range(0, n_samples, sample_chunk_size):
            s_end = min(s_start + sample_chunk_size, n_samples)
            sample_slice = slice(s_start, s_end)

            regions_dict = ag.get_many_regions(
                regions=batch_genes,
                sample_slice=sample_slice,
            )

            # Preallocate burden arrays for this gene batch
            sum_b = np.full((s_end - s_start, len(batch_genes)), np.nan, dtype=np.float32)
            max_b = np.full_like(sum_b, np.nan)
            top2_b = np.full_like(sum_b, np.nan)

            for g_idx, gene in enumerate(batch_genes):
                region = regions_dict.get(gene)
                if region is None:
                    continue

                # This returns (samples, 1)
                s, m, t = get_gene_burdens_numba(
                    region["genotypes"],
                    region["annotations"],
                    annotation_list=[anno],
                    max_burden=True,
                )

                sum_b[:, g_idx] = s[:, 0]
                max_b[:, g_idx] = m[:, 0]
                top2_b[:, g_idx] = t[:, 0]

                del region
                gc.collect()

            # Write to Zarr
            sum_burdens[sample_slice, gene_slice, anno_idx] = sum_b
            max_burdens[sample_slice, gene_slice, anno_idx] = max_b
            top2_burdens[sample_slice, gene_slice, anno_idx] = top2_b

            # Free memory
            del sum_b, max_b, top2_b
            gc.collect()
            
    print(f"Stored burdens in Zarr at {output_zarr}")
        
