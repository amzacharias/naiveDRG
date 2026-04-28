#!/bin/bash
#SBATCH --job-name=st18
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=st18.out
#SBTACH --error=st18.err

# Title: StringTie to assemble transcripts, pass 1
# Author: Amanda Zacharias
# Date: 2023-07-12
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load StdEnv/2020 stringtie/2.1.5

# Set variables
INPUT=absolutePath/mouNaiveSNI/2_align/aligned/18.sort.mrkdup.bam
OUT_GTF=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/pass1gtfs/18.gtf
REF_GTF=absolutePath/mouNaiveSNI/0_resources/gencode/gencode.vM24.primary_assembly.annotation.gtf

echo Stringtie Started at $(date +%T)
stringtie $INPUT -p 5 -G $REF_GTF -o $OUT_GTF

echo Job ended at $(date +%T)
