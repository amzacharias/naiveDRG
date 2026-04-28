#!/bin/bash
#SBATCH --job-name=run_seurat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --cpus-per-task=5
#SBATCH --mem=30MB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=run_seurat.out
#SBATCH --error=run_seurat.err
# Title: Run Seurat analysis
# Author: Amanda Zacharias
# Date: 2025-09-02
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
#
# Code -------------------------------------------
echo Job started at "$(date +%T)"
# Options

# Dependencies
module load StdEnv/2023 r/4.5.0
# Variables

# Body
echo ------------------------------- dea markers -------------------------------
Rscript dea_markers.R
echo ------------------------------- dea sex -------------------------------
Rscript dea_sex.R

echo Job ended at "$(date +%T)"