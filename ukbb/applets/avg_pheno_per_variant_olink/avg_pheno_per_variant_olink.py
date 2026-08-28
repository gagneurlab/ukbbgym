#!/usr/bin/env python3

"""
TO BUILD: dx build applets/avg_pheno_per_variant_olink/ --destination project-REDACTED:/applets/ -f
"""

import glob
import math
import os
import dxpy
import logging
import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

WORK_DIR = "PATH_TO_FILE"
IN_DIR   = f"{WORK_DIR}/in"
TMP_DIR  = f"{WORK_DIR}/out/tmp"
OUT_DIR  = f"{WORK_DIR}/out"
OUTPUT_NAME = "average_pheno_per_variant_olink.parquet"



def resolve_input(dxlink, local_name):
    local_path = os.path.join(IN_DIR, local_name)
    logger.info(f"Downloading {local_name}...")
    dxpy.download_dxfile(dxlink, local_path)
    return local_path


def load_associations(path):
    df = (
        pl.read_parquet(path)
        .rename({"gene": "gene_id"})
        .select(["gene_id"])
        .with_columns(phenotype=pl.col("gene_id") + "_olink")
    )
    logger.info(f"Associations loaded: {df.shape}")
    return df


def load_phenotypes(path, samples, olink_regions):
    """Load PROTRIDER-corrected Olink phenotypes (no INT). Intersects columns with whitelist regions."""
    phenos = pl.read_parquet(path)
    genes_to_use = list(set(phenos.columns).intersection(set(olink_regions)))
    phenos = phenos.select(["sample"] + genes_to_use)

    if samples is not None:
        phenos = phenos.filter(pl.col("sample").is_in(samples))

    long_phenos = (
        phenos
        .unpivot(
            index="sample",
            on=genes_to_use,
            variable_name="phenotype",
            value_name="pheno_value",
        )
        .with_columns(phenotype=pl.col("phenotype") + "_olink")
        .drop_nulls()
    )
    logger.info(f"Olink phenotypes loaded: {long_phenos.shape}")
    return long_phenos


def filter_variants(ann_path, genes, mac_threshold):
    filters = [pl.col("region").is_in(genes)]
    if mac_threshold is not None:
        filters.append(pl.col("mac_ukb") <= mac_threshold)
        logger.info(f"Filtering by gene_id and MAC <= {mac_threshold}")
    else:
        logger.info(f"Filtering by gene_id (no MAC filtering)")

    id_list = (
        pl.scan_parquet(ann_path)
        .filter(filters)
        .select("id")
        .unique()
        .collect()
    )
    logger.info(f"Variant id list: {id_list.shape}")
    return id_list


def build_carrier_genotypes(gt_path, id_list, samples):
    filters = [pl.col("gt") == 1]
    if samples is not None:
        filters.append(pl.col("sample").is_in(samples))

    long_gt = (
        pl.scan_parquet(gt_path)
        .select(["id", "sample", "gt"])
        .filter(filters)
        .join(id_list.lazy(), on="id", how="semi")
        .collect()
    )
    logger.info(f"Carrier genotypes loaded: {long_gt.shape}")
    return long_gt


def process_chunk(chunk_phenos, long_phenos, long_gt, tmp_path):
    (
        long_phenos.lazy()
        .filter(pl.col("phenotype").is_in(chunk_phenos))
        .join(long_gt.lazy(), on="sample", how="inner")
        .group_by(["id", "phenotype"])
        .agg(
            n_individuals=pl.len().cast(pl.Int32),
            mean_pheno_value=pl.col("pheno_value").mean().cast(pl.Float32),
            std_pheno_value=pl.col("pheno_value").std().cast(pl.Float32),
        )
        .with_columns(
            mean_pheno_value_rank=pl.col("mean_pheno_value")
                .rank(descending=True, method="max")
                .over("phenotype")
                .cast(pl.Float32),
        )
        .sink_parquet(tmp_path)
    )


def concat_and_filter(tmp_glob, ann_path, whitelist_df, out_path):
    """
    Memory-safe concat: filter each raw chunk to whitelist gene-trait pairs before combining.
    Avoids the ~92 GB intermediate that would result from concatenating all raw chunks first.
    """
    raw_files = sorted(glob.glob(tmp_glob))
    logger.info(f"Filtering and collecting {len(raw_files)} chunk files...")

    id_gene = (
        pl.scan_parquet(ann_path)
        .filter(pl.col("region").is_in(whitelist_df["gene_id"]))
        .select(["id", "region"])
        .rename({"region": "gene_id"})
        .unique()
        .collect()
    )
    logger.info(f"id_gene lookup: {id_gene.shape}")

    filtered_chunks = []
    for f in raw_files:
        chunk = (
            pl.scan_parquet(f)
            .join(id_gene.lazy(), on="id", how="inner")
            .join(whitelist_df.lazy(), on=["gene_id", "phenotype"], how="inner")
            .drop("gene_id")
            .unique(subset=["id", "phenotype"])
            .collect(engine="streaming")
        )
        filtered_chunks.append(chunk)
        logger.info(f"  {os.path.basename(f)}: {chunk.shape}")

    result = pl.concat(filtered_chunks)
    logger.info(f"Combined result: {result.shape}")
    result.write_parquet(out_path)
    logger.info(f"Final output written to {out_path}")


@dxpy.entry_point("main")
def main(
    associations_parquet,
    phenotypes_parquet,
    annotations_parquet,
    genotypes_parquet,
    pheno_chunk_size=100,
    mac_threshold=None,
    samples_txt=None,
):
    os.makedirs(IN_DIR, exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    assoc_path  = resolve_input(associations_parquet, "associations.parquet")
    pheno_path  = resolve_input(phenotypes_parquet,   "phenotypes.parquet")
    ann_path    = resolve_input(annotations_parquet,  "annotations.parquet")
    gt_path     = resolve_input(genotypes_parquet,    "gt_long.parquet")

    # Load associations (use all, no p-value filtering)
    whitelist_df = load_associations(assoc_path)

    # Load optional sample list
    samples = None
    if samples_txt is not None:
        sample_path = resolve_input(samples_txt, "samples.txt")
        with open(sample_path) as f:
            samples = [line.strip() for line in f if line.strip()]
        logger.info(f"Samples loaded from file: {len(samples)}")
    else:
        logger.info("No sample list provided; using all samples from phenotypes file")

    # Load Olink phenotypes (no INT)
    genes = whitelist_df["gene_id"].unique().to_list()
    long_phenos = load_phenotypes(pheno_path, samples, genes)

    # Filter variants by gene_id and MAC
    id_list = filter_variants(ann_path, genes, mac_threshold)

    # Build carrier genotype table (use all samples if none specified)
    gt_samples = samples if samples is not None else None
    long_gt = build_carrier_genotypes(gt_path, id_list, gt_samples)

    # Process phenotype chunks
    all_phenos = long_phenos["phenotype"].unique().to_list()
    num_phenos = len(all_phenos)
    logger.info(f"Processing {num_phenos} phenotypes in chunks of {pheno_chunk_size}...")

    for i in range(0, num_phenos, pheno_chunk_size):
        chunk_phenos = all_phenos[i: i + pheno_chunk_size]
        tmp_path = os.path.join(TMP_DIR, f"chunk_{i:04d}.parquet")
        logger.info(f"Chunk {i // pheno_chunk_size + 1}/{math.ceil(num_phenos / pheno_chunk_size)}: {len(chunk_phenos)} phenotypes")
        process_chunk(chunk_phenos, long_phenos, long_gt, tmp_path)

    # Concatenate with memory-safe pre-filtering
    out_path = os.path.join(OUT_DIR, OUTPUT_NAME)
    concat_and_filter(os.path.join(TMP_DIR, "*.parquet"), ann_path, whitelist_df, out_path)

    # Upload output
    uploaded = dxpy.upload_local_file(out_path)
    logger.info(f"Uploaded output: {uploaded.get_id()}")
    return {"appv_parquet": dxpy.dxlink(uploaded)}


dxpy.run()
