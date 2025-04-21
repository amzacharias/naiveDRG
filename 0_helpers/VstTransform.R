#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: VST Transformation
# Author: Amanda Zacharias
# Date: 2023-11-30
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(DESeq2) # 1.26.0


# Function -----------------------------------------
VstTransform <- function(counts, coldat) {
  #' Perform variance stabilizing transformation
  #' @param counts Gene counts matrix
  #' @param coldata Dataframe with metadata about samples
  #' @return Returns a transformed gene count matrix
  dds <- DESeqDataSetFromMatrix(round(counts), coldat, ~sampleGrps)
  dds <- estimateSizeFactors(dds)
  vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  return(vst)
}





