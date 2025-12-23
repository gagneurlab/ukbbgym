import polars as pl
import pyranges as pr
import numpy as np
from bgen import BgenWriter
import logging
import sys
import os

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)



def write_sample_file(sample_ids, output_path):
    """Writes the .sample file required by BGEN format."""
    logger.info(f"\nWriting .sample file to: {output_path}")
    header = "ID_1 ID_2 missing\n"
    second_header = "0 0 0\n"
    with open(output_path, 'w') as f:
        f.write(header)
        f.write(second_header)
        for sample_id in sample_ids:
            f.write(f"{sample_id} {sample_id} 0\n")
    logger.info("... .sample file written successfully.")

def write_bgen_from_scores(scores_df, gene_info, bgen_path):
    """
    Writes a BGEN file directly from a Polars DataFrame of scores
    using the correct bgen-py API keywords (varid, rsid, chrom, pos).
    """
    logger.info(f"\nWriting .bgen file to: {bgen_path}")
    
    sample_ids = scores_df.get_column("IID").to_list()
    n_samples = len(sample_ids)

    logger.info(f"{n_samples} samples found in the scores DataFrame.")
    
    # The BGEN writer works as a context manager
    with BgenWriter(bgen_path, n_samples=n_samples) as writer:
        # Iterate through each gene (column) in the scores dataframe
        for gene_name in scores_df.columns:
            if gene_name == "IID":
                continue
            logger.info(f"Processing gene: {gene_name}")  
            # Get the metadata for the current gene
            gene_meta = gene_info.filter(pl.col("gene_name") == gene_name)

            # Get the vector of scores for this gene
            scores = scores_df.get_column(gene_name).to_numpy()
            
            # --- Convert scores to probabilities (same logic as before) ---
            min_score, max_score = scores.min(), scores.max()
            if max_score > min_score:
                scaled_scores = (scores - min_score) / (max_score - min_score)
            else:
                scaled_scores = np.zeros_like(scores)

            probabilities = np.zeros((n_samples, 3))
            #probabilities[:, 0] = 1 - scaled_scores
            #probabilities[:, 1] = scaled_scores
            
            probabilities[:, 0] = scaled_scores/2
            probabilities[:, 1] = 0
            probabilities[:, 2] = 1 - (scaled_scores/2)


            if gene_meta.is_empty():
                logger.info(f"Warning: No metadata found for gene {gene_name}. Skipping it.")
                
            else:
                writer.add_variant(
                    varid = gene_name,                                  # Argument 1: varid
                    rsid = gene_name,                                   # Argument 2: rsid
                    chrom = str(gene_meta.get_column("Chromosome")[0]),      # Argument 3: chromosome
                    pos = gene_meta.get_column("Start")[0],               # Argument 4: position
                    alleles=[
                            "A",
                            "C",
                        ],                                               # Argument 6: allele2
                    genotypes = probabilities,
                    ploidy=2,
                    bit_depth=16,
                )

    print("... .bgen file written successfully.")


def generate_gene_metadata(gtf_file):
    """
    Creates a gene metadata DataFrame based on a GTF file.
    """
    gene_pos = pr.read_gtf(gtf_file)
    gene_pos = gene_pos[
        (gene_pos.Feature == "gene") & (gene_pos.gene_type == "protein_coding")
    ][["Chromosome", "Start", "End", "gene_id"]].as_df()

    gene_meta = pl.from_pandas(gene_pos)
    gene_meta = gene_meta.with_columns(
        pl.col('gene_id').str.split('.').list.first().alias('gene_name')
    )
    
    logger.info("... Gene metadata created successfully.")
    return gene_meta

# ==============================================================================
# 3. RUN THE PURE PYTHON FUNCTIONS
# ==============================================================================

OUT_FOLDER = "PATH_TO_FILE"
scores_df = pl.read_parquet("PATH_TO_FILE")
scores_df = scores_df.rename({'sample': 'IID'}).fill_null(0)
gtf_file = 'PATH_TO_FILE'

# Generate the gene_info_df automatically
gene_info_df = generate_gene_metadata(gtf_file)

sample_ids = scores_df.get_column('IID').to_list()
logger.info(f"{len(sample_ids)} samples found in the scores DataFrame.")


write_sample_file(sample_ids, OUT_FOLDER + "scores_bgen_lofteeHC.sample")
write_bgen_from_scores(scores_df, gene_info_df, OUT_FOLDER + "scores_bgen_lofteeHC.bgen")

logger.info("\nBGEN and SAMPLE files have been created directly from Python!")
logger.info("You can now run REGENIE Step 2.")