#!/bin/bash
#SBATCH --job-name=
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=
#SBTACH --error=

# Title: Run edgeR
# Author: Amanda Zacharias
# Date: 2023-08-11
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

# Variables
RSCRIPTNAME=
BASENAME=
PROJNAME=
PREFIX=
THRESH=
CANDPATH=
CWD=

# Execute R script
Rscript ../../${RSCRIPTNAME} \
--basename $BASENAME \
--projectName $PROJNAME \
--prefix $PREFIX \
--threshold $THRESH \
--candPath $CANDPATH \
--countPath ${CWD}/${PROJNAME}/2_dataPrep/${PREFIX}/cleanData/rawCounts.csv \
--workingDir $CWD
  
echo Job ended at $(date +'%T')
