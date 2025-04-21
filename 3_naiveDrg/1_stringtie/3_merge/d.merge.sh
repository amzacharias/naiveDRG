#!/bin/bash
#SBATCH --job-name=d.merge
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=10
#SBATCH --mem=3GB  # Job memory request
#SBATCH --time=0-2:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=d.merge.out
#SBTACH --error=d.merge.err

# Title: This script uses StringTie to merge gtfs
# Author: Amanda Zacharias
# Date: 2023-08-02
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo "job started"

# Load dependencies
module load StdEnv/2020 stringtie/2.1.5
   
# Set variables
BASEPATH=absolutePath/mouNaiveSNI
GTFS_LIST=${BASEPATH}/3_naiveDrg/1_stringtie/gtfLists/d.pass1List.txt
OUTPUT=${BASEPATH}/3_naiveDrg/1_stringtie/3_merge/d.merged.gtf
REF_GTF=${BASEPATH}/0_resources/gencode/gencode.vM24.primary_assembly.annotation.gtf

stringtie --merge -p 10 -o $OUTPUT -G $REF_GTF $GTFS_LIST

echo "job ended"
