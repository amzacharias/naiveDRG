#!/bin/bash
#SBATCH --job-name=gcd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=gcd.out
#SBTACH --error=gcd.err

# Title: Gff Compare
# Author: Amanda Zacharias
# Date: 2023-07-20
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load nixpkgs/16.09
module load gffcompare/0.11.6

REFGTFPATH=absolutePath/mouNaiveSNI/0_resources/gencode/gencode.vM24.primary_assembly.annotation.gtf
MERGEGTFPATH=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/3_merge/d.merged.gtf
OUTPREFIX=absolutePath/mouNaiveSNI/3_naiveDrg/1_stringtie/4_gffCompare/d

# Q
gffcompare -Q -T -R -r $REFGTFPATH -o $OUTPREFIX/Q $MERGEGTFPATH
# No Q
gffcompare -T -R -r $REFGTFPATH -o $OUTPREFIX/NoQ $MERGEGTFPATH

echo Job ended at $(date +%T)
