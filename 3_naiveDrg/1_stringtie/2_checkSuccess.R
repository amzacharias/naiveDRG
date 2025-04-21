#-------------------------------------------------
# Title: Check success of stringtie scripts
# Author: Amanda Zacharias
# Date: 2023-01-01
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# If job not yet completed successfully, write "sbatch <scriptn>" to a text file
# Run these jobs manually, the output file will just
# save you some typing and make it easier to monitor progress !
#
# Packages -----------------------------------------
library(dplyr) # 1.0.7

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "3_naiveDrg", "1_stringtie")
pass1gtfsDir <- file.path(baseDir, "pass1gtfs")
pass2gtfsDir <- file.path(baseDir, "pass2gtfs")
outFilePath <- file.path(baseDir, "jobsToRun.sh")

# Load data -----------------------------------------
# Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Make list of expected out files -----------------------------------------
makeExpectedList <- function(directory, coldata) {
  paths <- c()
  for (sample in coldata$sampleNum) { # loop through samples
    path <- file.path(
      directory,
      paste(sample, "gtf", sep = ".")
    )
    paths <- c(paths, path)
  } # finish looping through samples
  return(paths)
}
pass1ExpectedFiles <- makeExpectedList(pass1gtfsDir, coldata)
pass2ExpectedFiles <- makeExpectedList(pass2gtfsDir, coldata)

# Check files -----------------------------------------
# Check whether file exists and is sufficient in size
# If file isn't good, save to a list of lines to be written
checkFiles <- function(paths, pass_num, file_end, tool) {
  badLines <- c(paste("Pass number is", pass_num))
  for (path in paths) { # loop through samples
    reqSize <- 100000000
    if (file.exists(path) == FALSE | file.size(path) < reqSize) {
      cat("\n", path, file.info(path)$size)
      sample <- gsub(file_end, "", basename(path))
      newLine <- paste(
        "sbatch",
        paste(tool,
          paste(sample, "sh", sep = "."),
          sep = ""
        )
      )
      badLines <- c(badLines, newLine)
    } # end if statement
  } # end loop through samples
  return(badLines)
}
pass1Lines <- checkFiles(pass1ExpectedFiles, 1, ".gtf", "st")
pass2Lines <- checkFiles(pass2ExpectedFiles, 2, ".gtf", "st")

# Write lines -----------------------------------------
fileConn <- file(outFilePath)
writeLines(c(pass1Lines, pass2Lines), fileConn)
close(fileConn)
