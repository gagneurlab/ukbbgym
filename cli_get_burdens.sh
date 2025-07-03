#!/bin/bash

set -e  # Exit immediately if a command fails

# === Config ===
BASE_DIR="PATH_TO_FILE"
CONFIG_PATH="ukbgym/config_wgs_absplice2.yaml"
ASSOCIATIONS_DF_PATH="data_dir/absplice2_assocs.parquet"
OUTPUT_DIR="absplice2_78_small_genes"
GENE_CHUNK_SIZE="78"
SAMPLE_CHUNK_SIZE="40000"

# === Run the Python script ===
python scripts/get_burdens_chunky.py \
    --config-path "$BASE_DIR/$CONFIG_PATH" \
    --associations-df-path "$BASE_DIR/$ASSOCIATIONS_DF_PATH" \
    --output-dir "$BASE_DIR/$OUTPUT_DIR" \
    --only-snps \
    --na-mask \
    --overwrite \
    --gene-chunk-size "$GENE_CHUNK_SIZE" \
    --sample-chunk-size "$SAMPLE_CHUNK_SIZE" \
    > "$BASE_DIR/stdout.log" 2> "$BASE_DIR/stderr.log"

# === Upload results to DNAnexus ===
TARGET_FOLDER="project-REDACTED:/processed_data/ukbgym/burdens/$OUTPUT_DIR"

# Create directory recursively if not exists
dx mkdir -p "$TARGET_FOLDER"

# Upload files after folder exists
dx upload -r "$BASE_DIR/$OUTPUT_DIR" --path "$TARGET_FOLDER"

