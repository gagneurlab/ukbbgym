#!/bin/bash

set -e  # Exit immediately if a command fails

# === Config ===
BASE_DIR="/home/dnanexus"
CONFIG_PATH="ukbgym/config_wgs_cadd.yaml"
ASSOCIATIONS_DF_PATH="small_gene_assocs.pq"
OUTPUT_DIR="55_small_genes_onlySNP_cadd_annotations"
GENE_CHUNK_SIZE="55"
SAMPLE_CHUNK_SIZE="40000"

# === Run the Python script ===
python scripts/get_burdens_chunky.py \
    --config-path "$BASE_DIR/$CONFIG_PATH" \
    --associations-df-path "$BASE_DIR/$ASSOCIATIONS_DF_PATH" \
    --output-dir "$BASE_DIR/$OUTPUT_DIR" \
    --only-snps \
    --overwrite \
    --gene-chunk-size "$GENE_CHUNK_SIZE" \
    --sample-chunk-size "$SAMPLE_CHUNK_SIZE" \
    > "$BASE_DIR/stdout.log" 2> "$BASE_DIR/stderr.log"

# === Upload results to DNAnexus ===
TARGET_FOLDER="project-Gyp4fvjJg0yFZ374KvP9bGFJ:/processed_data/ukbgym/burdens/$OUTPUT_DIR"

# Create directory recursively if not exists
dx mkdir -p "$TARGET_FOLDER"

# Upload files after folder exists
dx upload -r "$BASE_DIR/$OUTPUT_DIR" --path "$TARGET_FOLDER"

