#!/bin/bash
#SBATCH --job-name=st11
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=st11.out
#SBTACH --error=st11.err

# Title: StringTie to assemble transcripts, pass 2
# Author: Amanda Zacharias
# Date: 2023-07-12
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load StdEnv/2020 stringtie/2.1.5

# Set variables
INPUT=absolutePath/mouNaiveSNI/2_align/aligned/11.sort.mrkdup.bam
OUT_GTF=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/pass2gtfs/11.gtf
REF_GTF=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/3_merge/d.merged.gtf
BALL=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/ballgown/d/11

echo Stringtie Started at $(date +%T)
stringtie $INPUT -b $BALL -e -p 5 -G $REF_GTF -o $OUT_GTF

echo Job ended at $(date +%T)
