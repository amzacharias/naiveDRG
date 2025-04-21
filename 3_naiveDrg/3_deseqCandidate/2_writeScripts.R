#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual bash scripts
# Author: Amanda Zacharias
# Date: 2023-08-11
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
projectName <- "3_naiveDrg"
prefix <- "gene"
rscriptNames <- list(
  "2vs14" = "0_2vs14.R"
)
thresholds <- list(
  "gene" = 0.015571
)

# Packages -----------------------------------------
library(dplyr) # 1.1.0


# Pathways -----------------------------------------
# Input ===========
baseDir <- file.path(getwd(), projectName, "3_deseqCandidate")
candidatesDir <- file.path(baseDir, "candidates", "cleanCandidates")
candPaths <- list.files(candidatesDir, full.names = TRUE)
names(candPaths) <- gsub(".csv", "", basename(candPaths))

# Output ===========
bashDir <- file.path(baseDir, prefix, "bash")
system(paste("mkdir", bashDir))

# Load data -----------------------------------------
baseScript <- as.character(
  read.table(file.path(baseDir, "2_baseScript.sh"),
    sep = "\n", blank.lines.skip = FALSE,
    comment.char = "", quote = "\'",
    stringsAsFactors = FALSE
  )$V1
)

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, rscriptPath, candPath, projName, prefix, thresh, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param rscriptPath Path to Rscript that will be run (string)
  #' @param candPath Path to candidate genes txt file (string)
  #' @param projName string; name of project folder
  #' @param prefix string; name of project sub-folder
  #' @param thresh double; threshold for non-specific filtering
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[25] <- paste(script[25], rscriptPath, sep = "")
  script[26] <- paste(script[26], name, sep = "")
  script[27] <- paste(script[27], projName, sep = "")
  script[28] <- paste(script[28], prefix, sep = "")
  script[29] <- paste(script[29], thresh, sep = "")
  script[30] <- paste(script[30], candPath, sep = "")
  script[31] <- paste(script[31], getwd(), sep = "")
  # Return
  return(script)
}
ProcessInfo <- function(candPaths, rscriptNames) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  # Loop through samples
  for (candName in names(candPaths)) {
    cat("\n", candName)
    for (rName in names(rscriptNames)) {
      cat("\n\t", rName)
      sampleLines <- ModifyScript(
        rscriptPath = rscriptNames[[rName]],
        candPath = candPaths[[candName]],
        name = paste(candName, rName, sep = "."),
        projName = projectName,
        prefix = prefix,
        script = baseScript,
        thresh = thresholds[[prefix]]
      )
      # Save
      outFilename <- paste(candName, rName, "sh", sep = ".")
      fileConn <- file(file.path(bashDir, outFilename))
      writeLines(sampleLines, fileConn)
      close(fileConn)
    } # end loop through rNames
  } # end loop through candNames
} # end function

# Execute -----------------------------------------
ProcessInfo(candPaths, rscriptNames)

# To Run
toRun <- tidyr::crossing(names(candPaths), names(rscriptNames))
toRunStr <- paste(toRun$`names(candPaths)`, toRun$`names(rscriptNames)`, "sh", sep = ".")
toRunLines <- paste("sbatch", toRunStr)

fileConn <- file(file.path(baseDir, prefix, "jobsToRun.sh"))
writeLines(toRunLines, fileConn)
close(fileConn)
