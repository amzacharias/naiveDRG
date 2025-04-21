#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write pass 2 scripts for stringtie
# Author: Amanda Zacharias
# Date: 2023-06-30
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "3_naiveDrg", "1_stringtie")
# Input ===========
coldataPath <- file.path("coldata.csv")
baseScriptPath <- file.path(baseDir, "5_baseScript2.sh")
mergeGtfPaths <- list.files(file.path(baseDir, "3_merge"),
  pattern = ".gtf", full.names = TRUE
)
names(mergeGtfPaths) <- gsub(".merged.gtf", "", basename(mergeGtfPaths))
alignedDir <- file.path(getwd(), "2_align", "aligned")
inputFilepaths <- list.files(alignedDir, full.names = TRUE)
inputFilepaths <- inputFilepaths[!grepl(".bai", inputFilepaths)]

# Output ===========
scriptsDir <- file.path(baseDir, "pass2indivScripts")
gtfsDir <- file.path(baseDir, "pass2gtfs")
ballDir <- file.path(baseDir, "ballgown")

system(paste("mkdir", scriptsDir, gtfsDir, ballDir))

# Load data -----------------------------------------
# Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Base script =======
baseScript <- read.table(
  baseScriptPath,
  sep = "\n", blank.lines.skip = FALSE,
  comment.char = "", quote = "\'", stringsAsFactors = FALSE
)$V1 %>%
  as.character()

# Function to get sample names -----------------------------------------
getNames <- function(paths) {
  names <- gsub(".sort.mrkdup.bam", "", basename(paths))
  return(names)
}

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, refPath, outPath, inPath, ballPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param refPath path to reference genome gtf file (string)
  #' @param outPath path to output gtf file (string)
  #' @param inPath path to input .bam file (string)
  #' @param ballPath path to output ballgown directory (string)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[25] <- paste(script[25], inPath, sep = "")
  script[26] <- paste(script[26], outPath, sep = "")
  script[27] <- paste(script[27], refPath, sep = "")
  script[28] <- paste(script[28], ballPath, sep = "")
  # Return
  return(script)
}
ProcessInfo <- function(coldat) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  for (set in unique(coldat$tissue)) {
    cat("\n", set)
    # Tissue-specific folders
    setBallDir <- file.path(ballDir, set)
    system(paste("mkdir", setBallDir))
    # Loop through tissue-specific samples
    setColdata <- coldat %>% subset(tissue == set)
    for (idx in 1:nrow(setColdata)) {
      # Path info
      sampleNum <- setColdata$sampleNum[idx]
      samplePathIdx <- match(sampleNum, getNames(inputFilepaths))
      samplePath <- inputFilepaths[samplePathIdx]
      # Modify base script
      sampleLines <- ModifyScript(
        name = sampleNum,
        refPath = mergeGtfPaths[[set]],
        outPath = file.path(gtfsDir, paste(sampleNum, "gtf", sep = ".")),
        inPath = samplePath,
        ballPath = file.path(setBallDir, sampleNum),
        script = baseScript
      )
      # Save
      filename <- paste("st", sampleNum, ".sh", sep = "")
      fileConn <- file(file.path(scriptsDir, filename))
      writeLines(sampleLines, fileConn)
      close(fileConn)
    } # end loop through samples in a tissue
  } # end loop through tissues
} # end function

# Execute ===================================
ProcessInfo(coldata)
