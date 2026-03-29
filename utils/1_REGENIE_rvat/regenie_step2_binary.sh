#!/bin/bash
#SBATCH --job-name=rvat_regenie  # Job name
#SBATCH --output=./logs/regenie-%j.stdout       # Output log file
#SBATCH --error=./logs/regenie-%j.stderr         # Error log file
#SBATCH --cpus-per-task=64          # Number of CPU cores per task
#SBATCH --mem=128G                   # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

# --- Command to run ---
# sbatch ./regenie_step2.sh

# Input Files
COVARIATE_FILE="/s/project/geno2pheno/BFuncRVP/data/regenie/regenie_input/covariates.txt"
BGEN_FILE="/s/project/deeprvat/ukb_gym/regenie/scores_bgen.bgen"
SAMPLE_FILE="/s/project/deeprvat/ukb_gym/regenie/sample_bgen.sample"
PHENOTYPE_FILE="/s/project/geno2pheno/BFuncRVP/data/regenie/regenie_input/jurgens_binary_phenos.txt"
STEP1_PRED_FILE="/s/project/deeprvat/ukb_gym/regenie/binary_traits/prs.list_pred.list"

echo "Reading phenotype names from header of: $PHENOTYPE_FILE"

# This awk command reads the header (NR==1), starts from the 3rd field (i=3)
# to skip FID and IID, and prints each subsequent field followed by a space.
# The 'for' loop in bash will iterate over these space-separated names.
ALL_TRAITS=$(awk 'NR==1 {for(i=3; i<=NF; i++) printf "%s ", $i} {exit}' "$PHENOTYPE_FILE")
echo "Found $(echo $ALL_TRAITS | wc -w) phenotypes to process."

# Output Prefix and Directories
OUTPUT_DIR="/s/project/deeprvat/ukb_gym/regenie/binary_traits/rvat"
OUTPUT_PREFIX="${OUTPUT_DIR}/step2_outputs"            # Base name for output files
# EXPECTED_OUTPUT_FILE="${OUTPUT_PREFIX}_prs.list" # For reference, matches Snakemake output definition

# Log Files
LOG_DIR="/s/project/deeprvat/ukb_gym/regenie/logs"
LOG_STDOUT="${LOG_DIR}/regenie_step2_binary.stdout"
LOG_STDERR="${LOG_DIR}/regenie_step2_binary.stderr"

# Parameters
THREADS=64                                      # Number of threads
REGENIE_STEP2_BSIZE=1000
# Corresponds to regenie_config_step1.get("options", [])
# Add any extra regenie command-line options here as separate elements
REGENIE_EXTRA_OPTIONS=(
    # "--example-option"
    # "--another-option value" # If an option takes a value, keep them together if possible, or handle separately
)

# Temporary Directory
TEMP_DIR="/scratch/tmp/regenie_rvat_tmp_binary/"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TEMP_DIR"

# --- Script Logic ---
echo "Starting REGENIE Step 2..."
echo "Input BGEN: $BGEN_FILE"
echo "Sample File: $SAMPLE_FILE"
echo "Phenotype File: $PHENOTYPE_FILE"
echo "Covariate File: $COVARIATE_FILE"
echo "Output Prefix: $OUTPUT_PREFIX"
echo "Log Directory: $LOG_DIR"
echo "Threads: $THREADS"
echo "Extra Options: ${REGENIE_EXTRA_OPTIONS[*]}"
echo "Temporary Directory: $TEMP_DIR"

echo "Running REGENIE..."


# --- Create a DYNAMIC output prefix for this specific trait ---
# This prevents each run from overwriting the previous one's results.
OUTPUT_PREFIX="${OUTPUT_DIR}/step2_results_${TRAIT_NAME}"

echo "Output files will be prefixed with: $OUTPUT_PREFIX"

# Execute the REGENIE command for the current trait
regenie \
    --step 2 \
    --bt \
    --bgen "$BGEN_FILE" \
    --sample "$SAMPLE_FILE" \
    --phenoFile "$PHENOTYPE_FILE" \
    --covarFile "$COVARIATE_FILE" \
    --pred "$STEP1_PRED_FILE" \
    --bsize "$REGENIE_STEP2_BSIZE" \
    --out "$OUTPUT_PREFIX" \
    > "$LOG_STDOUT" 2> "$LOG_STDERR"

echo "REGENIE step 2 finished successfully. Cleaning up temporary directory..."
rm -rf "$TEMP_DIR"

echo "Script completed."
