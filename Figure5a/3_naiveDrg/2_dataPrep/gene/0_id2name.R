#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make id2name file
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

# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "3_naiveDrg")
# Input ===========
t2gPath <- file.path(baseDir, "1_stringtie", "6_isoformAnalyzeR", "d.t2g.csv")

# Output ===========
id2namePath <- file.path(baseDir, "2_dataPrep", prefix, "id2name.csv")

# Load data -----------------------------------------
t2g <- read.csv(t2gPath, stringsAsFactors = FALSE, row.names = 1) %>%
  dplyr::select(-c("isoform_id")) %>%
  distinct()

# Save data -----------------------------------------
write.csv(t2g, id2namePath)
