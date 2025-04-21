#!/bin/bash
#SBATCH --job-name=indexGenHisat2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB  # Job memory request
#SBATCH --time=0-10:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=indexGenHisat2.out
#SBTACH --error=indexGenHisat2.err

# Building HGFM index with transcript
# Amanda Zacharias, July 7th, 2022

# Reference link used to make this script: 
# http://daehwankimlab.github.io/hisat2/howto/ 

echo "job started"

# Load HISAT2 and dependencies
module load StdEnv/2020
module load hisat2/2.2.1

# Assign variables
BASEPATH=absolutePath/mouNaiveSNI
RESOURCEDIR=${BASEPATH}/0_resources/gencode
IDXDIR=${BASEPATH}/2_align/indexGencode
mkdir $IDXDIR

# Reference genome paths
FNA_REF=${RESOURCEDIR}/GRCm38.primary_assembly.genome.fa
GENOMENAME=`basename $FNA_REF .fa`
GTF_REF=${RESOURCEDIR}/gencode.vM24.primary_assembly.annotation.gtf

# Output index paths
INDEX_NAME=${IDXDIR}/idx
SPLICESITE_PATH=${IDXDIR}/${GENOMENAME}.ss
EXON_PATH=${IDXDIR}/${GENOMENAME}.exon

# Extract splice sites and exons
echo "Started extracting splice sites"
hisat2_extract_splice_sites.py $GTF_REF > $SPLICESITE_PATH
echo "Started extracting exons"
hisat2_extract_exons.py $GTF_REF > $EXON_PATH

# Build index
echo "Started building the index"
hisat2-build -p 10 -f $FNA_REF --ss $SPLICESITE_PATH --exon $EXON_PATH $INDEX_NAME

echo "job ended"
