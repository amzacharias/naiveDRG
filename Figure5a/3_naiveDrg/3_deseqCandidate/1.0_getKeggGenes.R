#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Download genes in KEGG pathways of interest
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
options(stringsAsFactors = FALSE)
prefix <- "3_naiveDrg"

# Packages -----------------------------------------
library(biomaRt) # 2.42.1


# Pathways -----------------------------------------
baseDir <- getwd()
edgeRCandDir <- file.path(baseDir, prefix, "3_deseqCandidate", "candidates", "rawCandidates")
keggCandDir <- file.path(baseDir, prefix, "3_deseqCandidate", "candidates", "keggRaw")
system(paste("mkdir", edgeRCandDir, keggCandDir))

# Load data -----------------------------------------
# Download pathway and genes info from KEGG ==========
keggMouseLink <- "http://rest.kegg.jp/link/mmu"
keggIDs <- list(
  "IonChannels" = "br:mmu04040"
)
for (pathID in names(keggIDs)) {
  system(paste("wget -P", keggCandDir, file.path(keggMouseLink, keggIDs[[pathID]])))
}

# Read in as dataframes ==========
pathDfs <- lapply(
  1:length(keggIDs),
  function(n) read.table(file.path(keggCandDir, keggIDs[[n]]), col.names = c("path", "gene"))
)
names(pathDfs) <- names(keggIDs)

# Get genes -----------------------------------------
queryGeneIDs <- lapply(
  1:length(pathDfs),
  function(n) data.frame(entrezgene_id = gsub(".*:", "", pathDfs[[n]]$gene))
)
names(queryGeneIDs) <- names(pathDfs)

# Convert ids to gene names -----------------------------------------
ensembl <- useEnsembl(biomart = "ensembl")
# listDatasets(ensembl)
mouMart <- useDataset(
  dataset = "mmusculus_gene_ensembl", mart = ensembl
) # Mouse genes (GRCm39)

RunBiomart <- function(geneIDs) {
  bmQueryRes <- getBM(
    filters = "entrezgene_id",
    values = geneIDs$entrezgene_id,
    attributes = c(
      "entrezgene_id", "ensembl_gene_id", "external_gene_name",
      "external_gene_source"
    ),
    mart = mouMart,
    useCache = FALSE
  ) # useCache = FALSE to get rid of error
  bmMergedRes <- merge(geneIDs, bmQueryRes, byå = "entrezgene_id")
  colnames(bmMergedRes) <- c(
    "entrezgene_id", "ensembl_gene_id", "external_gene_name",
    "external_gene_source"
  )
  return(bmMergedRes)
}

bmResList <- lapply(
  1:length(queryGeneIDs),
  function(n) RunBiomart(queryGeneIDs[[n]])
)
names(bmResList) <- names(queryGeneIDs)

# Save ------------------------------
for (pathName in names(bmResList)) {
  toWrite <- bmResList[[pathName]] %>%
    dplyr::rename("gene_name" = "external_gene_name")
  write.csv(
    toWrite,
    file.path(edgeRCandDir, paste(keggIDs[[pathName]], "csv", sep = "."))
  )
}
