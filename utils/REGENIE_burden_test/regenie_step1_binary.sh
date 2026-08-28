#!/bin/bash
#SBATCH --job-name=regenie  # Job name
#SBATCH --output=./logs/regenie-%j.stdout       # Output log file
#SBATCH --error=./logs/regenie-%j.stderr         # Error log file
#SBATCH --cpus-per-task=64         # Number of CPU cores per task
#SBATCH --mem=256G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

# --- Command to run ---
# sbatch -p urgent ./regenie_step1.sh

# Input Files
BGEN_FILE="PATH_TO_FILE"
SNP_LIST="PATH_TO_FILE"
COVARIATE_FILE="PATH_TO_FILE"
PHENOTYPE_FILE="PATH_TO_FILE"
# PHENOTYPE_FILE="PATH_TO_FILE"

# For Samples
SAMPLE_LIST="PATH_TO_FILE"
echo "Detected Sample subset file: $SAMPLE_LIST"

# For Covariates:
COVAR_LIST=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s%s", $i, (i<NF ? "," : "")} {exit}' "$COVARIATE_FILE")
echo "Detected Covariates: $COVAR_LIST"

# For Phenotypes:
PHENO_LIST=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s%s", $i, (i<NF ? "," : "")} {exit}' "$PHENOTYPE_FILE")
# Only print the first few detected phenotypes to keep the log clean
echo "Detected Phenotypes (first few): $(echo $PHENO_LIST | cut -d, -f1-5)..."

# Output Prefix and Directories
# OUTPUT_DIR="PATH_TO_FILE"
OUTPUT_DIR="PATH_TO_FILE"
OUTPUT_PREFIX="${OUTPUT_DIR}/prs.list"            # Base name for output files
# EXPECTED_OUTPUT_FILE="${OUTPUT_PREFIX}_prs.list" # For reference, matches Snakemake output definition

# Log Files
LOG_DIR="PATH_TO_FILE"
LOG_STDOUT="${LOG_DIR}/regenie_step1_binary.stdout"
LOG_STDERR="${LOG_DIR}/regenie_step1_binary.stderr"

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
TEMP_DIR="/scratch/tmp/regenie_prs_binary_tmp/"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TEMP_DIR"

# --- Script Logic ---
echo "Starting REGENIE Step 1..."
echo "Input BGEN: $BGEN_FILE"
echo "Input SNP List: $SNP_LIST"
echo "Phenotype File: $PHENOTYPE_FILE"
echo "Sample subset file: $SAMPLE_LIST"
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
    --bt \
    --print-prs	\
    --bgen "$BGEN_FILE" \
    --extract "$SNP_LIST" \
    --phenoFile "$PHENOTYPE_FILE" \
    --covarFile "$COVARIATE_FILE" \
    --phenoColList "$PHENO_LIST" \
    --keep "$SAMPLE_LIST" \
    --covarColList "$COVAR_LIST" \
    --bsize "$REGENIE_STEP1_BSIZE" \
    --threads "$THREADS" \
    --lowmem \
    --lowmem-prefix "${TEMP_DIR}/regenie_prs_binary" \
    "${REGENIE_EXTRA_OPTIONS[@]}" \
    --out "$OUTPUT_PREFIX" \
    > "$LOG_STDOUT" 2> "$LOG_STDERR"

# If REGENIE completes successfully (due to '&&' implicit in set -e and explicit chaining if used), clean up
echo "REGENIE step 1 finished successfully. Cleaning up temporary directory..."
rm -rf "$TEMP_DIR"

echo "Script completed."

# Optional: Check if the expected output file exists
# if [ -f "$EXPECTED_OUTPUT_FILE" ]; then
#     echo "Expected output file found: $EXPECTED_OUTPUT_FILE"
# else
#     echo "Warning: Expected output file NOT found: $EXPECTED_OUTPUT_FILE" >&2 # Write warning to stderr
# fi
