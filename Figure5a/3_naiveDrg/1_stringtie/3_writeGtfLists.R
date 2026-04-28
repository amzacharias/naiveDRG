#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write list of paths to GTFs
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
# Input ===========
coldataPath <- file.path("coldata.csv")
baseDir <- file.path(getwd(), "3_naiveDrg", "1_stringtie")
gtfs1Dir <- file.path(baseDir, "pass1gtfs")
gtfs2Dir <- file.path(baseDir, "pass2gtfs")

# Output ===========
gtfListsDir <- file.path(baseDir, "gtfLists")
system(paste("mkdir", gtfListsDir))

# Load data -----------------------------------------
# Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Write paths -----------------------------------------
writePaths <- function(coldat, gtfDir, newFilename) {
  lines <- c()
  for (idx in 1:nrow(coldat)) {
    sample <- coldat$sampleNum[idx]
    path <- file.path(gtfDir, paste(sample, "gtf", sep = "."))
    if (grepl("2", newFilename)) {
      lines <- c(lines, paste(sample, path))
    } else {
      lines <- c(lines, path)
    }
  } # end loop through samples
  cat("writing lines with", length(lines), "lines\n")
  write(lines, file.path(gtfListsDir, newFilename))
}
writePaths(coldata %>% subset(tissue == "d"), gtfs1Dir, "d.pass1List.txt")
writePaths(coldata %>% subset(tissue == "d"), gtfs2Dir, "d.pass2List.txt")
