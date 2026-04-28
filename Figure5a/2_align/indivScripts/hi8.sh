#!/bin/bash
#SBATCH --job-name=hi8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=5
#SBATCH --mem=30GB  # Job memory request #SBATCH --tmp=20GB  # Temporary memory
#SBATCH --time=1-10:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=hi8.out
#SBTACH --error=hi8.err

# Aligning reads with hisat2.
# Amanda Zacharias
# July 7th, 2022

echo Job started at $(date +%T)

# Load dependencies
module load StdEnv/2020 samtools/1.10 hisat2/2.2.1 

# Set variables
PAIRED_END_1=absolutePath/mouNaiveSNI/0_data/F20FTSUSA1339_MOUpfnaT/CleanData/8/8_1.fq.gz
PAIRED_END_2=absolutePath/mouNaiveSNI/0_data/F20FTSUSA1339_MOUpfnaT/CleanData/8/8_2.fq.gz
INDEX=absolutePath/mouNaiveSNI/2_align/indexGencode/idx
SAM_PATH=absolutePath/mouNaiveSNI/2_align/aligned/8.sam
OUT_PATH=absolutePath/mouNaiveSNI/2_align/aligned/8
SUMMARY_PATH=absolutePath/mouNaiveSNI/2_align/summaries/8.txt

# Start alignment
echo ALignment Started at $(date +%T)
hisat2 -p 10 -x $INDEX -1 $PAIRED_END_1 -2 $PAIRED_END_2 \
  --dta --sensitive --no-discordant --no-mixed  \
  --summary-file $SUMMARY_PATH --time --verbose -S $SAM_PATH

# Process with samtools
echo Samtools processing started at $(date +%T)
samtools view -b -@ 10 $SAM_PATH > ${OUT_PATH}.bam
rm $SAM_PATH
echo collate started at $(date +%T)
samtools collate -@ 10 -o ${OUT_PATH}.col.bam  ${OUT_PATH}.bam ${OUT_PATH}_tmpcol
rm ${OUT_PATH}.bam
echo fixmate started at $(date +%T)
samtools fixmate -m -@ 10 ${OUT_PATH}.col.bam  ${OUT_PATH}.fix.bam 
rm ${OUT_PATH}.col.bam
echo sort started at $(date +%T)
samtools sort -@ 10 -T ${OUT_PATH}_sort -o ${OUT_PATH}.sort.bam ${OUT_PATH}.fix.bam 
rm ${OUT_PATH}.fix.bam
echo markdup started at $(date +%T)
samtools markdup -@ 10 -T ${OUT_PATH}_tmpmrk -s ${OUT_PATH}.sort.bam ${OUT_PATH}.sort.mrkdup.bam 
rm ${OUT_PATH}.sort.bam
echo index started at $(date +%T)
samtools index -b -@ 10 ${OUT_PATH}.sort.mrkdup.bam  ${OUT_PATH}.sort.mrkdup.bam.bai

echo Job ended at $(date +%T)
