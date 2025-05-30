#!/bin/bash
#SBATCH --job-name=fq4_1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB  # Job memory request
#SBATCH --time=0-12:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=fq4_1.out
#SBTACH --error=fq4_1.err

# Title: This script uses fastqc to check qualtiy of sequencing reads
# Author: Amanda Zacharias
# Date: 2023-07-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 

echo Job started at $(date +%T)

module load StdEnv/2020
module load nixpkgs/16.09
module load fastqc/0.11.9
module list

INDATAPATH=absolutePath/mouNaiveSNI/0_data/F20FTSUSA1339_MOUpfnaT/CleanData/4/4_1.fq.gz
OUTDIR=absolutePath/mouNaiveSNI/1_qcSeqReads/fastqcOut
fastqc -f fastq -o $OUTDIR $INDATAPATH

echo Job ended at $(date +%T)
