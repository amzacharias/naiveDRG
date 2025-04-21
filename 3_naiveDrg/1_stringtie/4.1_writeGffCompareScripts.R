#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write GFF Compare Scripts
# Author: Amanda Zacharias
# Date: 2023-07-20
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
refGtfPath <- file.path(
  getwd(), "0_resources", "gencode",
  "gencode.vM24.primary_assembly.annotation.gtf"
)
baseScriptPath <- file.path(baseDir, "5_baseGffCompare.sh")
mergeGtfsDir <- file.path(baseDir, "3_merge")
mergeGtfs <- list.files(mergeGtfsDir, pattern = ".gtf", full.names = TRUE)
names(mergeGtfs) <- gsub(".merged.gtf", "", basename(mergeGtfs))

baseScriptPath <- file.path(baseDir, "4_baseGffCompare.sh")

# Output ===========
gffCompareDir <- file.path(baseDir, "4_gffCompare")
system(paste(
  "mkdir", gffCompareDir,
  file.path(gffCompareDir, "d")
))

# Load data -----------------------------------------
# Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

baseScript <- read.table(
  baseScriptPath,
  sep = "\n", blank.lines.skip = FALSE,
  comment.char = "", quote = "\'",
  stringsAsFactors = FALSE
)$V1 %>%
  as.character()

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, refPath, outPath, inPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param refPath path to reference genome gtf file (string)
  #' @param outPath path to output files, including filename prefix (string)
  #' @param inPath path to input gtf file (string)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[25] <- paste(script[25], refPath, sep = "")
  script[26] <- paste(script[26], inPath, sep = "")
  script[27] <- paste(script[27], outPath, sep = "")
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
    # Modify base script
    sampleLines <- ModifyScript(
      name = set,
      refPath = refGtfPath,
      outPath = file.path(gffCompareDir, set),
      inPath = mergeGtfs[[set]],
      script = baseScript
    )
    # Save
    filename <- paste("gc", set, ".sh", sep = "")
    fileConn <- file(file.path(gffCompareDir, filename))
    writeLines(sampleLines, fileConn)
    close(fileConn)
  } # end loop through tissues
}

# Execute ===================================
ProcessInfo(coldata)
