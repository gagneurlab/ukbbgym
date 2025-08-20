#!/bin/bash
#SBATCH --job-name=regenie  # Job name
#SBATCH --partition=noninterruptive
#SBATCH --output=./logs/regenie-%j.stdout       # Output log file
#SBATCH --error=./logs/regenie-%j.stderr         # Error log file
#SBATCH --cpus-per-task=32         # Number of CPU cores per task
#SBATCH --mem=256G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

CHUNK="$1"

# Input Files
BGEN_FILE="PATH_TO_FILE"
SNP_LIST="PATH_TO_FILE"
COVARIATE_FILE="PATH_TO_FILE"
PHENOTYPE_FILE="PATH_TO_FILE"


# For Covariates:
COVAR_LIST=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s%s", $i, (i<NF ? "," : "")} {exit}' "$COVARIATE_FILE")
echo "Detected Covariates: $COVAR_LIST"

# For Phenotypes:
PHENO_LIST=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s%s", $i, (i<NF ? "," : "")} {exit}' "$PHENOTYPE_FILE")
# Only print the first few detected phenotypes to keep the log clean
echo "Detected Phenotypes (first few): $(echo $PHENO_LIST | cut -d, -f1-5)..."


# Output Prefix and Directories
OUTPUT_DIR="PATH_TO_FILE"
OUTPUT_PREFIX="${OUTPUT_DIR}/output_step1_chunk_${CHUNK}"            # Base name for output files
# EXPECTED_OUTPUT_FILE="${OUTPUT_PREFIX}_prs.list" # For reference, matches Snakemake output definition

# Log Files
LOG_DIR="PATH_TO_FILE"
LOG_STDOUT="${LOG_DIR}/regenie_step1.stdout"
LOG_STDERR="${LOG_DIR}/regenie_step1.stderr"

# Parameters
THREADS=32                                      # Number of threads
REGENIE_STEP1_BSIZE=1000                        # Corresponds to regenie_step1_bsize (Example value!)
# Corresponds to regenie_config_step1.get("options", [])
# Add any extra regenie command-line options here as separate elements
REGENIE_EXTRA_OPTIONS=(
    # "--example-option"
    # "--another-option value" # If an option takes a value, keep them together if possible, or handle separately
)

# Temporary Directory
TEMP_DIR="/scratch/tmp/regenie_prs_tmp_${CHUNK}/"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TEMP_DIR"

# --- Script Logic ---
echo "Starting REGENIE Step 1..."
echo "Input BGEN: $BGEN_FILE"
echo "Input SNP List: $SNP_LIST"
echo "Phenotype File: $PHENOTYPE_FILE"
echo "Covariate File: $COVARIATE_FILE"
echo "Output Prefix: $OUTPUT_PREFIX"
echo "Log Directory: $LOG_DIR"
echo "Threads: $THREADS"
echo "Block Size: $REGENIE_STEP1_BSIZE"
echo "Extra Options: ${REGENIE_EXTRA_OPTIONS[*]}"
echo "Temporary Directory: $TEMP_DIR"

echo "Running REGENIE..."

# Execute the REGENIE command with redirection
# Note: Ensure the regenie executable is in your PATH
# REMOVE phenoColList and covarColList if you don't want to save prediction
regenie \
    --step 1 \
    --print-prs	\
    --bgen "$BGEN_FILE" \
    --extract "$SNP_LIST" \
    --phenoFile "$PHENOTYPE_FILE" \
    --covarFile "$COVARIATE_FILE" \
    --phenoColList "$PHENO_LIST" \
    --covarColList "$COVAR_LIST" \
    --bsize "$REGENIE_STEP1_BSIZE" \
    --threads "$THREADS" \
    --lowmem \
    --lowmem-prefix "${TEMP_DIR}/regenie_prs" \
    "${REGENIE_EXTRA_OPTIONS[@]}" \
    --out "$OUTPUT_PREFIX" \
    > "$LOG_STDOUT" 2> "$LOG_STDERR"

# If REGENIE completes successfully (due to '&&' implicit in set -e and explicit chaining if used), clean up
echo "REGENIE step 1 finished successfully. Cleaning up temporary directory..."
rm -rf "$TEMP_DIR"

echo "Script completed."




# Input Files
BGEN_FILE="PATH_TO_FILE"
SAMPLE_FILE="PATH_TO_FILE"
STEP1_PRED_FILE="PATH_TO_FILE"

echo "Reading phenotype names from header of: $PHENOTYPE_FILE"

# This awk command reads the header (NR==1), starts from the 3rd field (i=3)
# to skip FID and IID, and prints each subsequent field followed by a space.
# The 'for' loop in bash will iterate over these space-separated names.
ALL_TRAITS=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s ", $i} {exit}' "$PHENOTYPE_FILE")

echo "Found $(echo $ALL_TRAITS | wc -w) phenotypes to process."



# COVARIATE_FILE="PATH_TO_FILE"
# PHENOTYPE_FILE="PATH_TO_FILE"     
# SAMPLE_FILE="/path/to/keep_samples.txt"       # Corresponds to input.sample_file (if uncommented)

# Output Prefix and Directories
OUTPUT_DIR="PATH_TO_FILE"

# Log Files
LOG_DIR="PATH_TO_FILE"
LOG_STDOUT="${LOG_DIR}/regenie_step2.stdout"
LOG_STDERR="${LOG_DIR}/regenie_step2.stderr"

# Parameters
THREADS=32                                      # Number of threads
REGENIE_STEP1_BSIZE=1000                        # Corresponds to regenie_step1_bsize (Example value!)
# Corresponds to regenie_config_step1.get("options", [])
# Add any extra regenie command-line options here as separate elements
REGENIE_EXTRA_OPTIONS=(
    # "--example-option"
    # "--another-option value" # If an option takes a value, keep them together if possible, or handle separately
)



# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"



for TRAIT_NAME in $ALL_TRAITS; do

    echo "============================================================"
    echo "         STARTING REGENIE STEP 2 FOR TRAIT: $TRAIT_NAME"
    echo "============================================================"

    # --- Create a DYNAMIC output prefix for this specific trait ---
    # This prevents each run from overwriting the previous one's results.
    OUTPUT_PREFIX="${OUTPUT_DIR}/step2_results_${TRAIT_NAME}"

    echo "Output files will be prefixed with: $OUTPUT_PREFIX"

    # Execute the REGENIE command for the current trait
    regenie \
      --step 2 \
      --bgen "$BGEN_FILE" \
      --sample "$SAMPLE_FILE" \
      --phenoFile "$PHENOTYPE_FILE" \
      --covarFile "$COVARIATE_FILE" \
      --phenoCol "$TRAIT_NAME" \
      --pred "$STEP1_PRED_FILE" \
      --bsize 200 \
      --out "$OUTPUT_PREFIX" \
        > "$LOG_STDOUT" 2> "$LOG_STDERR"

    echo "--- Finished processing trait: $TRAIT_NAME ---"
    echo

done

echo "============================================================"
echo "All phenotype processing complete."
echo "============================================================"

