#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Outlier removal
# Author: Amanda Zacharias
# Date: 2023-08-02
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
prefix <- "gene"
options(ggrepel.max.overlaps = Inf, box.padding = 0)
source("0_helpers/OutlierDetect.R")
source("0_helpers/ExtractOutliers.R")
source("0_helpers/VstTransform.R")
source("0_helpers/MakePCA.R")
source("0_helpers/MakeHeatmap.R")

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(Biobase) # ‘2.46.0’
library(arrayQualityMetrics) # 3.42.0
library(DESeq2) # 1.26.0
library(pheatmap) # 1.0.12
# devtools::install_github("marcosci/cividis")
library(cividis) # 0.2.0
library(ggplot2) # 3.3.6
library(ggrepel) # 0.9.1
library(stringr) # 1.4.0

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "3_naiveDrg")
# Input ===========
coldataPath <- file.path("coldata.csv")
countsPath <- file.path(baseDir, "1_stringtie", "6_isoformAnalyzeR", "d.geneCounts.csv")

# Output ===========
dataPrepDir <- file.path(baseDir, "2_dataPrep", prefix)
qcB4OutDir <- file.path(dataPrepDir, "qcB4Out")
qcAftOutDir <- file.path(dataPrepDir, "qcAftOut")
cleanDataDir <- file.path(dataPrepDir, "cleanData")
plotsDir <- file.path(dataPrepDir, "plots")
heatmapsDir <- file.path(plotsDir, "heatmaps")
pcaDir <- file.path(plotsDir, "pca")

system(paste(
  "mkdir", qcB4OutDir, qcAftOutDir,
  cleanDataDir, plotsDir, heatmapsDir,
  pcaDir
))

# Load data -----------------------------------------
coldata <- read.csv(coldataPath, row.names = 1) %>%
  mutate(sampleGrps = as.factor(sampleGrps))
rownames(coldata) <- coldata$sampleNum
counts <- read.csv(countsPath, row.names = 1, check.names = FALSE) %>%
  tibble::column_to_rownames("gene_id") %>%
  dplyr::select(-c("gene_name"))

# Match order of samples -----------------------------------------
countsMatchDf <- counts %>% dplyr::select(as.character(coldata$sampleNum))

# Matrix -----------------------------------------
countsMat <- countsMatchDf %>% data.matrix()

# Raw QC -----------------------------------------
rawOutliers <- OutlierDetect(
  counts = countsMat,
  coldat = coldata,
  outDirname = file.path(qcB4OutDir, "raw"),
  intGroup = "ztTime"
)
saveRDS(rawOutliers, file = file.path(dataPrepDir, "rawOutliers.rds"))

# Normalize counts -----------------------------------------
Normalize <- function(data, coldat) {
  #' Perform TMM normalization of gene counts
  #'
  #' @param data A dataframe / matrix of gene counts
  #' @param coldat Metadata about samples
  #' @return Returns a normalized dataframe
  #' @example
  #' normDf <- Normalize(counts, FALSE)
  #' normLogDf <- Normalize(counts, TRUE)
  dds <- DESeqDataSetFromMatrix(round(data.matrix(data)), coldat, ~sampleGrps)
  ddsSE <- estimateSizeFactors(dds)
  return(counts(ddsSE, normalized = TRUE))
}
normCountsMat <- Normalize(countsMat, coldata)

# Norm QC -----------------------------------------
normOutliers <- OutlierDetect(
  counts = normCountsMat,
  coldat = coldata,
  outDirname = file.path(qcB4OutDir, "norm"),
  intGroup = "sampleGrps"
)
saveRDS(normOutliers, file = file.path(dataPrepDir, "normOutliers.rds"))

# Remove outliers round 1 -----------------------------------------
names(rawOutliers) <- paste("raw", names(rawOutliers)) # only run this once!
names(normOutliers) <- paste("norm", names(normOutliers)) # only run this once!
round1Outliers <- append(rawOutliers, normOutliers)
round1ToRemove <- ExtractOutliers(round1Outliers)

capture.output(round1ToRemove, file = file.path(dataPrepDir, "round1ToRemove.txt"))
saveRDS(round1ToRemove, file = file.path(dataPrepDir, "round1ToRemove.rds"))
round1ToRemove <- readRDS(file.path(dataPrepDir, "round1ToRemove.rds"))

# No outliers!

# VST transformation -----------------------------------------
vstRaw <- VstTransform(counts, coldata)

# Get pca plot data -----------------------------------------
rawPcaDf <- GetPcaData(
  vstDat = vstRaw, 
  removeVec = NA, 
  toPlotVars = c("sampleGrps", "ztTime")
  )

# Sample PCA plot -----------------------------------------
MakePCA(
  pcaDat = rawPcaDf,
  newFilename = "raw.pdf",
  newPath = pcaDir
)

# Sample heatmap plot -----------------------------------------
MakeSampleHeatmap(
  countMat = assay(vstRaw),
  varOfInterest = "ztTime",
  coldat = coldata, 
  newFilename = "raw.pdf",
  newPath = heatmapsDir
)

# Save files -----------------------------------------
write.csv(counts,
  file = file.path(cleanDataDir, "rawCounts.csv")
)
write.csv(normCountsMat,
  file = file.path(cleanDataDir, "normCounts.csv")
)
write.csv(coldata,
  file = file.path(cleanDataDir, "coldata.csv")
)
write.csv(assay(vstRaw),
  file = file.path(cleanDataDir, "vstCounts.csv")
)

# Save image -------------------------------------------
save.image(file = file.path(dataPrepDir, "outlierRemoval.RData"))

