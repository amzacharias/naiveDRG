#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Non-specific filtering
# Author: Amanda Zacharias
# Date: 2023-08-07
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
prefix <- "gene"
source("0_helpers/filtering.R")

# Packages -----------------------------------------
library(dplyr) # 1.1.0


# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "3_naiveDrg", "2_dataPrep", prefix)
# Input ===========
countsPaths <- list.files(file.path(baseDir, "cleanData"),
  full.names = TRUE, pattern = "Counts.csv"
)
names(countsPaths) <- gsub("Counts.csv", "", basename(countsPaths))

# Output ===========
filtData <- file.path(baseDir, "filtData")
plotsDir <- file.path(baseDir, "plots")
system(paste("mkdir", filtData))

# Load data -----------------------------------------
countsList <- lapply(
  1:length(countsPaths),
  function(n) read.csv(countsPaths[[n]], row.names = 1, check.names = FALSE)
)
names(countsList) <- names(countsPaths)

# Nonspecific filtering -----------------------------------------
# Plot =======
GetMadCutoff(countsList$norm, plotsDir)

# Remove lowly expressed features =======
# Using variance quantiles
mads <- GetMads(countsList$norm)
quantile(mads$mad, probs = seq(0.1, 0.3, 0.05))
# 10%      15%      20%      25%      30%
# 0.000000 0.000000 0.000000 0.015571 2.189552
filtThreshold <- 0.015571
toKeep <- rownames(subset(mads, mad >= filtThreshold))

filtList <- lapply(
  1:length(countsList),
  function(n) {
    data.frame(countsList[[n]][rownames(countsList[[n]]) %in% toKeep, ],
      check.names = FALSE
    )
  }
)
names(filtList) <- names(countsList)

# What are the differences? ======
cat(
  cat("Before filtering: ", nrow(countsList$norm), "\n"),
  cat("After filtering: ", nrow(filtList$norm), "\n")
)
# Before filtering:  38942
# After filtering:  29206

# Save files -----------------------------------------
for (dfName in names(filtList)) {
  write.csv(
    filtList[[dfName]],
    file.path(filtData, paste(dfName, "Counts.csv", sep = ""))
  )
}

# Save image ------------------------------------------------------------------
save.image(file = file.path(baseDir, "filtering.RData"))
