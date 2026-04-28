#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: edgeR
# Author: Amanda Zacharias
# Date: 2023-08-08
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
options(ggrepel.max.overlaps = Inf, box.padding = 0)

# Packages -----------------------------------------
library(optparse) # 1.7.3
library(dplyr) # 1.1.0
library(tibble) # 3.1.8
library(DESeq2) # 1.26.0
library(ggplot2) # 3.4.1
library(ggrepel) # 0.9.3

# Run optparser -----------------------------------------
# Define arguments ==================
optionList <- list(
  make_option(c("-w", "--workingDir"),
    type = "character", default = getwd(),
    help = "the project working directory path", metavar = "character"
  ),
  make_option(c("-p", "--projectName"),
    type = "character", default = NA,
    help = "name of folder for analysis", metavar = "character"
  ),
  make_option(c("-x", "--prefix"),
    type = "character", default = NA,
    help = "name of sub-folder for analysis", metavar = "character"
  ),
  make_option(c("-c", "--countPath"),
    type = "character", default = NA,
    help = "path to your count matrix", metavar = "character"
  ),
  make_option(c("-n", "--basename"),
    type = "character", default = NA,
    help = "the base filename for this analysis' outputs", metavar = "character"
  ),
  make_option(c("-t", "--threshold"),
    type = "double", default = NA,
    help = "non-specific filtering MAD threshold", metavar = "double"
  ),
  make_option(c("-a", "--candPath"),
    type = "character", default = NA,
    help = "path to your candidates list", metavar = "character"
  )
)
# Get parameters ==================
optParser <- OptionParser(option_list = optionList)
opt <- parse_args(optParser)

# There is no parameter checking!!!!!
cat(
  "\nBase filename:", opt$basename,
  "\nProject name:", opt$projectName,
  "\nPrefix:", opt$prefix,
  "\nThreshold for filtering:", opt$threshold,
  "\n"
)

# Set the working directory -----------------------------------------
setwd(opt$workingDir)

# Helper functions -----------------------------------------
source("0_helpers/filtering.R")
source("0_helpers/MakeVolcano.R")

# Pathways -----------------------------------------
# Input ===========
coldataPath <- file.path(opt$projectName, "2_dataPrep", opt$prefix, "cleanData", "coldata.csv")
id2namePath <- file.path(opt$projectName, "2_dataPrep", opt$prefix, "id2name.csv")

# Output ===========
baseDir <- file.path(opt$projectName, "3_deseqCandidate", opt$prefix)
dfsDir <- file.path(baseDir, "dataframes")
baseDfsDir <- file.path(dfsDir, opt$basename)
rDataDir <- file.path(baseDir, "rData")
plotsDir <- file.path(baseDir, "plots")
basePlotsDir <- file.path(plotsDir, opt$basename)
system(paste(
  "mkdir", dfsDir, baseDfsDir, rDataDir,
  plotsDir, basePlotsDir
))

# Load data -----------------------------------------
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = FALSE)
id2name <- read.csv(id2namePath, row.names = 1, stringsAsFactors = FALSE)
counts <- read.csv(opt$countPath,
  row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE
)
candidates <- read.csv(opt$candPath, row.names = 1, stringsAsFactors = FALSE)
candidatesId2Name <- candidates %>%
  left_join(id2name, multiple = "all") %>%
  distinct()

# Execute DeSeq2 ########################################
# Normalize -----------------------------------------
cat("\n", rep("-", 10), "Normalizing", rep("-", 10), "\n")
dds <- DESeqDataSetFromMatrix(round(data.matrix(counts)), coldata, ~sampleGrps)
ddsSE <- estimateSizeFactors(dds)

# Dispersion estimation -----------------------------------------
cat("\n", rep("-", 10), "Dispersion Estimation", rep("-", 10), "\n")
ddsDisp <- estimateDispersions(ddsSE)

pdf(file.path(basePlotsDir, "dispersion.pdf"),
  width = 6, height = 6
)
plotDispEsts(ddsDisp,
  ylab = "Dispersion", xlab = "Mean of normalized counts",
  main = opt$basename
)
dev.off()

# Filtering -----------------------------------------
cat("\n", rep("-", 10), "Filtering", rep("-", 10), "\n")
GetMadCutoff(counts(ddsDisp, normalized = TRUE), basePlotsDir)
mads <- GetMads(counts(ddsDisp, normalized = TRUE))
toKeep <- rownames(subset(mads, mad >= opt$threshold))

ddsFilt <- ddsDisp[toKeep, ]
nrow(ddsFilt)

# Extract only candidates -----------------------------------------
cat("\n", rep("-", 10), "Extracting Candidates", rep("-", 10), "\n")
idCol <- "isoform_id"
if (opt$prefix == "gene") {
  idCol <- "gene_id"
}

ddsCand <- ddsFilt[match(candidatesId2Name[[idCol]], rownames(ddsFilt), nomatch = 0), ]

candGenesIncluded <- unique(candidatesId2Name$gene_name[
  match(rownames(ddsCand), candidatesId2Name[[idCol]], nomatch = 0)
])
candGenesNotIncluded <- unique(
  candidatesId2Name$gene_name[
    match(
      setdiff(candidatesId2Name[[idCol]], rownames(ddsCand)),
      candidatesId2Name[[idCol]],
      nomatch = 0
    )
  ]
)
cat(
  "\nNumber of candidates:", nrow(candidates),
  "\nNumber of features after subset:", length(rownames(ddsCand)),
  "\nNumber of unique features after subset:", length(unique(rownames(ddsCand))),
  "\nCandidate genes included:", candGenesIncluded,
  "\nCandidate genes not included:", candGenesNotIncluded,
  "\nCandidate genes not included at all:", setdiff(candGenesNotIncluded, candGenesIncluded),
  "\n"
)

# Differential expression analysis -----------------------------------------
cat("\n", rep("-", 10), "DE Test", rep("-", 10), "\n")
ddsWald <- nbinomWaldTest(ddsCand)
cont <- c("sampleGrps", "naivezt14d", "naivezt2d")
coefs <- "sampleGrps_naivezt14d_vs_naivezt2d"
res <- results(
  ddsWald,
  contrast = cont,
  independentFilter = FALSE,
  pAdjustMethod = "bonferroni",
  alpha = 0.05
)
summary(res)

# Format results -----------------------------------------
cat("\n", rep("-", 10), "Formatting Results", rep("-", 10), "\n")
resDf <- data.frame(res@listData, row.names = res@rownames)
if (opt$prefix == "transcript") {
  resDf <- resDf %>%
    arrange(padj) %>%
    tibble::rownames_to_column("isoform_id") %>%
    left_join(id2name, by = "isoform_id")
} else {
  resDf <- resDf %>%
    arrange(padj) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(id2name, by = "gene_id")
}
sigResDf <- resDf %>% subset(padj < 0.05)
cat("\nNumber of Sig:", nrow(sigResDf), "Out of:", nrow(resDf))

# Save results -----------------------------------------
cat("\n", rep("-", 10), "Saving Dataframes", rep("-", 10), "\n")
write.csv(resDf, file.path(baseDfsDir, "resDf.csv"))
write.csv(sigResDf, file.path(baseDfsDir, "sigResDf.csv"))

# Plot results -----------------------------------------
cat("\n", rep("-", 10), "Plotting", rep("-", 10), "\n")
# Volcano plot ========
MakeVolcano(
  df = resDf,
  showName = TRUE,
  newTitle = opt$basename,
  newPath = basePlotsDir,
  idColumn = idCol
)
MakeVolcano(
  df = resDf,
  showName = FALSE,
  newTitle = opt$basename,
  newPath = basePlotsDir,
  idColumn = idCol
)

# Save image -----------------------------------------
cat("\n", rep("-", 10), "Saving RData", rep("-", 10), "\n")
save.image(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))
