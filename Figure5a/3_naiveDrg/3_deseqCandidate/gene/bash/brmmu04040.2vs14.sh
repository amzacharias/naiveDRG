#!/bin/bash
#SBATCH --job-name=brmmu04040.2vs14
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=brmmu04040.2vs14.out
#SBTACH --error=brmmu04040.2vs14.err

# Title: Run edgeR
# Author: Amanda Zacharias
# Date: 2023-08-11
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

# Variables
RSCRIPTNAME=0_2vs14.R
BASENAME=brmmu04040.2vs14
PROJNAME=3_naiveDrg
PREFIX=gene
THRESH=0.015571
CANDPATH=absolutePath/mouNaiveSNI/3_naiveDrg/3_deseqCandidate/candidates/cleanCandidates/brmmu04040.csv
CWD=absolutePath/mouNaiveSNI

# Execute R script
Rscript ../../${RSCRIPTNAME} \
--basename $BASENAME \
--projectName $PROJNAME \
--prefix $PREFIX \
--threshold $THRESH \
--candPath $CANDPATH \
--countPath ${CWD}/${PROJNAME}/2_dataPrep/${PREFIX}/cleanData/rawCounts.csv \
--workingDir $CWD
  
echo Job ended at $(date +%T)
