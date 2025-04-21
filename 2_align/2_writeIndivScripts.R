#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual Hisat2 scripts
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
baseDir <- file.path(getwd(), "2_align")
# Input ========
baseScriptPath <- file.path(baseDir, "2_baseScript.sh")
inputFilepaths <- list.files(
  file.path(getwd(), "0_data", "F20FTSUSA1339_MOUpfnaT", "CleanData"),
  recursive = TRUE, full.names = TRUE, pattern = "fq.gz"
) %>%
  as.vector() # vector for `match()` to work
refGenIdxPath <- file.path(baseDir, "indexGencode", "idx")

# Output =======
scriptsDir <- file.path(baseDir, "indivScripts")
alignedDir <- file.path(baseDir, "aligned")
summariesDir <- file.path(baseDir, "summaries")
system(paste("mkdir", scriptsDir, alignedDir, summariesDir))

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
  names <- gsub("_[0-9].fq.gz", "", basename(paths))
  return(names)
}

# Write scripts -----------------------------------------
for (idx in 1:nrow(coldata)) {
  sample <- coldata$sampleNum[idx]
  cat(sample, "\n")
  # Assume that rev path immediately follows the fwd path in vector
  samplePathIdxs <- c(
    match(sample, getNames(inputFilepaths)),
    match(sample, getNames(inputFilepaths)) + 1
  )
  samplePaths <- inputFilepaths[samplePathIdxs]
  fwdPath <- samplePaths[grep("_1", samplePaths)]
  revPath <- samplePaths[grep("_2", samplePaths)]

  # Modify base script
  sampleLines <- baseScript # reset list of lines for every sample
  sampleLines[2] <- paste(sampleLines[2], sample, sep = "")
  sampleLines[9] <- paste(sampleLines[9], paste(sample, "out", sep = "."), sep = "")
  sampleLines[10] <- paste(sampleLines[10], paste(sample, "err", sep = "."), sep = "")
  sampleLines[22] <- paste(sampleLines[22], fwdPath, sep = "")
  sampleLines[23] <- paste(sampleLines[23], revPath, sep = "")
  sampleLines[24] <- paste(sampleLines[24], file.path(refGenIdxPath), sep = "")
  sampleLines[25] <- paste(sampleLines[25], file.path(
    alignedDir, paste(sample, "sam", sep = ".")
  ), sep = "")
  sampleLines[26] <- paste(sampleLines[26], file.path(alignedDir, sample), sep = "")
  sampleLines[27] <- paste(sampleLines[27], file.path(
    summariesDir, paste(sample, "txt", sep = ".")
  ), sep = "")

  # Write the file
  # https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r
  filename <- paste("hi", sample, sep = "")
  filename_sh <- paste(filename, "sh", sep = ".")
  fileConn <- file(file.path(scriptsDir, filename_sh))
  writeLines(sampleLines, fileConn)
  close(fileConn)
} # end idx loop (looping through samples)
