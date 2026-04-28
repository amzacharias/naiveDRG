#!/bin/bash
#SBATCH --job-name=gc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=gc
#SBTACH --error=gc

# Title: Gff Compare
# Author: Amanda Zacharias
# Date: 2023-07-20
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load nixpkgs/16.09
module load gffcompare/0.11.6

REFGTFPATH=
MERGEGTFPATH=
OUTPREFIX=

# Q
gffcompare -Q -T -R -r $REFGTFPATH -o $OUTPREFIX/Q $MERGEGTFPATH
# No Q
gffcompare -T -R -r $REFGTFPATH -o $OUTPREFIX/NoQ $MERGEGTFPATH

echo Job ended at $(date +'%T')
