#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Use IsoformAnalyzeR to convert transcript abundance --> gene counts
# Author: Amanda Zacharias
# Date: 2023-07-06
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load StdEnv/2020 r/4.2.1
# View section "Rescue StringTie Annotation and Extract Gene Count Matrix" in vignette
# https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#importing-data-into-r
# The newer the version of IsoformSwitchAnalyzeR, the better !
# To run in the terminal: 2 options
# R_LIBS_USER=absolutePath/R/x86_64-redhat-linux-gnu-library/4.2.1 R
# OR
# R_LIBS_USER=absolutePath/R/x86_64-redhat-linux-gnu-library/4.2.1 Rscript isoformAnalyzeR.R
#
# Options -----------------------------------------

# Packages -----------------------------------------
library(dplyr) # 1.1.0
# BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR) # 1.17.4

# Pathways -----------------------------------------
# Warning: using an absolute path!
setwd("absolutePath/mouNaiveSNI")
stringtieDir <- file.path("3_naiveDrg", "1_stringtie")
# Input ===========
coldataPath <- file.path("coldata.csv")
gtfPaths <- list.files(file.path(stringtieDir, "3_merge"), full.names = TRUE, pattern = ".gtf")
names(gtfPaths) <- gsub(".merged.gtf", "", basename(gtfPaths))
inputDirs <- file.path(stringtieDir, "ballgown", names(gtfPaths))
names(inputDirs) <- basename(inputDirs)

# Output ===========
isoDir <- file.path(stringtieDir, "6_isoformAnalyzeR")

# Load data -----------------------------------------
# Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)
coldata <- coldata %>% subset(tissue == "d" & naiveVsSni == "naive")
coldataList <- list("d" = coldata %>% subset(tissue == "d"))

# Format for isoformSwitchAnalyzeR
isoColdataList <- lapply(
  1:length(coldataList),
  function(n) {
    data.frame(
      sampleID = coldataList[[n]]$sampleNum,
      condition = coldataList[[n]]$ztTime
    )
  }
)
names(isoColdataList) <- names(coldataList)

# Get counts function -----------------------------------------
GetCounts <- function(inputDir, readlength, coldat, gtfPath) {
  #' Use IsoformSwitchAnalyzeR functions to get gene counts from transcript abundances
  #'
  #' @param inputDir A string path to the folder with ballgown files from StringTie.
  #' @param readlength The numeric length of sequencing reads.
  #' @param coldat Metadata dataframe explaining what each sample is.
  #' @param gtfPath A string path to the reference gtf used for assembly.
  #' @returns A named list.
  #' @examples
  #' countObj <- GetCounts('./ballgown', 150, df, './3_merge/merged.gtf')
  quant <- importIsoformExpression(
    parentDir = inputDir,
    readLength = readlength,
    addIsofomIdAsColumn = TRUE
  )
  switchQuant <- importRdata(
    isoformCountMatrix = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix = coldat,
    isoformExonAnnoation = gtfPath,
    showProgress = TRUE
  )
  geneCounts <- extractGeneExpression(
    switchQuant,
    extractCounts = TRUE # set to FALSE for abundances
  )
  return(list("switchQuant" = switchQuant, "geneCounts" = geneCounts))
}

# Execute function -----------------------------------------
countObjs <- list()
for (set in names(isoColdataList)) {
  cat("\n", set, "\n")
  # Get counts
  countObjs[[set]] <- GetCounts(
    inputDir = inputDirs[[set]],
    readlength = 150,
    coldat = isoColdataList[[set]],
    gtfPath = gtfPaths[[set]]
  )
  # Save files
  write.csv(
    countObjs[[set]]$geneCounts,
    file.path(isoDir, paste(set, "geneCounts.csv", sep = "."))
  )
  saveRDS(countObjs[[set]]$switchQuant,
    file = file.path(isoDir, paste(set, "switchQuant.rds", sep = "."))
  )
}

# Transcript id to gene id dataframe -----------------------------------------
GetT2G <- function(switchQuant) {
  #' Extract isoform_id to gene_id/gene_name mapping from switchQuant object
  #'
  #' @param switchQuant A switchQuant object from isoformSwitchAnalyzeR
  #' @returns A dataframe.
  #' @note Need a separate t2g for each dataframe because, for example,
  #' MSTRG.1 != MSTRG.1 in a different dataset.
  #' @examples
  #' t2g <- GetT2G(countObjs[[1]]$switchQuant)
  t2g <- switchQuant$isoformFeatures %>%
    dplyr::select(c("isoform_id", "gene_id", "gene_name")) %>%
    distinct(isoform_id, gene_id, gene_name)
  return(t2g)
}
t2gs <- lapply(
  1:length(countObjs),
  function(n) GetT2G(countObjs[[n]]$switchQuant)
)
names(t2gs) <- names(countObjs)

# Save
for (set in names(t2gs)) {
  write.csv(t2gs[[set]], file.path(isoDir, paste(set, "t2g.csv", sep = ".")))
}

# Output messages -----------------------------------------
# d
# Step 1 of 3: Identifying which algorithm was used...
# The quantification algorithm used was: StringTie
# Found 12 quantification file(s) of interest
# Step 2 of 3: Reading data...
# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12
# Step 3 of 3: Normalizing abundance values (not counts) via edgeR...
# Done
#
# Step 1 of 7: Checking data...
# Please note that some condition names were changed due to names not suited for modeling in R.
# Step 2 of 7: Obtaining annotation...
# importing GTF (this may take a while)...
# 29853 ( 19%) isoforms were removed since they were not expressed in any samples.
# Step 3 of 7: Fixing StringTie gene annoation problems...
# 10807 isoforms were assigned the ref_gene_id and gene_name of their associated gene_id.
# This was only done when the parent gene_id were associated with a single ref_gene_id/gene_name.
# 3251 isoforms were assigned the ref_gene_id and gene_name of the most similar
# annotated isoform (defined via overlap in genomic exon coordinates).
# This was only done if the overlap met the requriements
# indicated by the three fixStringTieViaOverlap* arguments.
# We were unable to assign 86 isoforms (located within annotated genes) to a known ref_gene_id/gene_name.
# These were removed to enable analysis of the rest of the isoform from within the merged genes.
# 1653 gene_ids which were associated with multiple ref_gene_id/gene_names
# were split into mutliple genes via their ref_gene_id/gene_names.
# 38622 genes_id were assigned their original gene_id instead of the StringTie gene_id.
# This was only done when it could be done unambiguous.
# Step 4 of 7: Calculating gene expression and isoform fractions...
# Step 5 of 7: Merging gene and isoform expression...
# |======================================================================| 100%
# Step 6 of 7: Making comparisons...
# |======================================================================| 100%
# Step 7 of 7: Making switchAnalyzeRlist object...
# The GUESSTIMATED number of genes with differential isoform usage are:
#   comparison estimated_genes_with_dtu
# 1 X14 vs X20                115 - 192
# 2  X14 vs X8                122 - 204
# 3  X20 vs X8                 95 - 159
# Done
#
# Warning messages:
#   1: In importRdata(isoformCountMatrix = quant$counts, isoformRepExpression = quant$abundance,  :
#                       We found 510 (0.32%) unstranded transcripts.
#                     These were removed as unstranded transcripts cannot be analysed
#                     2: In importRdata(isoformCountMatrix = quant$counts, isoformRepExpression = quant$abundance,  :
#                                         No CDS annotation was found in the GTF files meaning ORFs could not be annotated.
#                                       (But ORFs can still be predicted with the analyzeORF() function)
