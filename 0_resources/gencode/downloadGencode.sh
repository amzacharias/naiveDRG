#!/bin/bash
#SBATCH --job-name=downloadGencode
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=150GB  # Job memory request
#SBATCH --time=00-6:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=downloadGencode.out
#SBTACH --error=downloadGencode.err

# Download reference from genocde
cd absolutePath/mouNaiveSNI/0_resources/gencode

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz

gunzip -v *.gtf.gz
gunzip -v *.fa.gz
