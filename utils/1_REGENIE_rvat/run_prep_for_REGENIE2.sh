#!/bin/bash
#SBATCH --job-name=regenie  # Job name
#SBATCH --output=./logs/regenie-%j.stdout       # Output log file
#SBATCH --error=./logs/regenie-%j.stderr         # Error log file
#SBATCH --cpus-per-task=32         # Number of CPU cores per task
#SBATCH --mem=256G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

python prep_files_for_REGENIE2.py
