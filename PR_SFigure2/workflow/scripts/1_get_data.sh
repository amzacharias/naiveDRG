#!/bin/bash
#SBATCH --job-name=get_data
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --cpus-per-task=1
#SBATCH --mem=20MB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=get_data.out
#SBATCH --error=get_data.err
# Title: Download scRNA data from Renthal lab
# Author: Amanda Zacharias
# Date: 2025-07-22
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
#
# Code -------------------------------------------
echo Job started at "$(date +%T)"
# Options

# Renthal Data ====================
RENTHALLINK="https://www.dropbox.com/scl/fo/x9c173p8p0wxze15c29ln/AIEbaG7XNg9HncmOubsPnJQ/DRG_neurons_release.Rds?rlkey=gsg6vcp5sm8e88s36ffml4uur&dl=1"
DATADIR="../../data"
wget --quiet -P $DATADIR "$RENTHALLINK"
mv ${DATADIR}/"DRG_neurons_release.Rds?rlkey=gsg6vcp5sm8e88s36ffml4uur&dl=1" ${DATADIR}/"DRG_neurons_release.Rds"

echo Job ended at "$(date +%T)"