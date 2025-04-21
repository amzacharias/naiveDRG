#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Prepare candidates
# Author: Amanda Zacharias
# Date: 2023-08-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
prefix <- "3_naiveDrg"

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(stringr) # 1.5.0
library(org.Mm.eg.db) # 3.10.0

# Pathways -----------------------------------------
candidatesDir <- file.path(prefix, "3_deseqCandidate", "candidates")
# Input ===========
rawCandidatesDir <- file.path(candidatesDir, "rawCandidates")

# Output ===========
cleanCandidatesDir <- file.path(candidatesDir, "cleanCandidates")
system(paste("mkdir", cleanCandidatesDir))

# Load data -----------------------------------------
# Paths ====
candidateFiles <- list.files(rawCandidatesDir, full.names = TRUE, pattern = ".")
names(candidateFiles) <- gsub(".csv", "", basename(candidateFiles))

# Load =====
dfsList <- lapply(
  1:length(candidateFiles),
  function(n) read.csv(candidateFiles[[n]], stringsAsFactors = FALSE, row.names = 1)
)
names(dfsList) <- names(candidateFiles)

# Unify column name for genes -----------------------------------------
# SYMBOL to gene_name
# Keep only columns of interest, helps minimize duplicates
for (dfName in names(dfsList)) {
  cat("\n", dfName)
  if (grepl("GO", dfName) == TRUE) {
    cat("\tGO")
    dfsList[[dfName]] <- dfsList[[dfName]] %>%
      dplyr::rename("gene_name" = "SYMBOL") %>%
      dplyr::select(c("GOALL", "gene_name", "GENENAME", "ENSEMBL"))
  }
}

# To sentence -----------------------------------------
# Also ensure there are no repeats
toSList <- lapply(
  1:length(dfsList),
  function(n) {
    dfsList[[n]] <- dfsList[[n]] %>%
      mutate(gene_name = str_to_sentence(gene_name)) %>%
      dplyr::select(gene_name) %>%
      unique()
  }
)
names(toSList) <- names(dfsList)

# Save -----------------------------------------
for (dfName in names(toSList)) {
  cleanDfName <- gsub(":", "", dfName)
  write.csv(
    toSList[[dfName]],
    file.path(cleanCandidatesDir, paste(cleanDfName, "csv", sep = "."))
  )
}
