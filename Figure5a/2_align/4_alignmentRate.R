#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Summarize alignment rates
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
library(stringr) # 1.4.0
library(ggplot2) # 3.3.6

# Pathways -----------------------------------------
# Input ===========
toolDir <- file.path(getwd(), "2_align")
summariesDir <- file.path(toolDir, "summaries")
logPaths <- list.files(summariesDir, full.names = TRUE)

# Output ===========
ratesDir <- file.path(toolDir, "rates")
system(paste("mkdir", ratesDir))

# Load data -----------------------------------------
## Coldata =======
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)
coldata.naiveDrg <- coldata %>%
    subset((tissue == "d") & (naiveVsSni == "naive"))

## Subset logPaths by only coldata -----------------------------------------
logPaths.naiveDrg <- logPaths[
  match(coldata.naiveDrg$sampleNum, as.numeric(sub(".txt", "", basename(logPaths))))
]

# Get rates -----------------------------------------
GetRates <- function(logPaths) {
  # Get rates
  filenames <- c()
  uniqList <- c()
  multiList <- c()
  overallList <- c()
  for (idx in 1:length(logPaths)) {
    filepath <- logPaths[idx]
    filename <- gsub(".txt", "", basename(filepath))
    filenames <- c(filenames, filename)
    # unique rate info on 4th line
    uniq <- system(paste("sed '4q;d'", filepath, sep = " "), intern = TRUE)
    uniqNums <- str_extract(string = uniq, pattern = "(?<=\\().*(?=\\))")
    uniqNums <- as.numeric(substr(uniqNums, 1, nchar(uniqNums) - 1)) # remove %
    uniqList <- c(uniqList, uniqNums)
    # overall mapping rate on 6th line
    overall <- system(paste("sed '6q;d'", filepath, sep = " "), intern = TRUE)
    overallNums <- as.numeric(substr(overall, 1, 5))
    overallList <- c(overallList, overallNums)
  } # finish looping through samples
  # add to list
  ratesDfs <- data.frame(
    filename = filenames,
    overall = overallList,
    uniq = uniqList
  )
  cat("Done!", "\n")
  return(ratesDfs)
} # finish function
ratesDf <- GetRates(logPaths)
ratesDf.naiveDrg <- GetRates(logPaths.naiveDrg)

# Get statistics -----------------------------------------
GetStats <- function(ratesDf) {
  overallMin <- min(ratesDf$overall)
  overallMean <- mean(ratesDf$overall)
  overallMax <- max(ratesDf$overall)
  overallMedian <- median(ratesDf$overall)

  uniqMin <- min(ratesDf$uniq)
  uniqMean <- mean(ratesDf$uniq)
  uniqMax <- max(ratesDf$uniq)
  uniqMedian <- median(ratesDf$uniq)

  statsDf <- data.frame(
    groups = c("overall", "unique"),
    min = c(overallMin, uniqMin),
    mean = c(overallMean, uniqMean),
    max = c(overallMax, uniqMax),
    median = c(overallMedian, uniqMedian)
  )
  return(statsDf)
}
statsDf <- GetStats(ratesDf)
statsDf.naiveDrg <- GetStats(ratesDf.naiveDrg)

# Write dataframes -----------------------------------------
write.csv(ratesDf, file.path(ratesDir, "hisatRates.csv"))
write.csv(statsDf, file.path(ratesDir, "hisatStats.csv"))

write.csv(ratesDf.naiveDrg, file.path(ratesDir, "hisatRates.naiveDrg.csv"))
write.csv(statsDf.naiveDrg, file.path(ratesDir, "hisatStats.naiveDrg.csv"))

# Visualize stats -----------------------------------------
# Boxplot =========
MakeBoxplot <- function(df, newFilename) {
  boxplot <- df %>%
    dplyr::select(-"filename") %>%
    stack() %>%
    dplyr::rename("value" = "values", "variable" = "ind") %>%
    ggplot(aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    geom_point(
      alpha = 0.5, shape = 21,
      show.legend = FALSE, position = position_jitterdodge()
    ) +
    scale_y_continuous(name = "Alignment rate", breaks = seq(0, 100, 10), limits = c(50, 100)) +
    scale_x_discrete(name = "Alignment rate type", labels = c("Overall", "Unique")) +
    scale_fill_viridis_d(
      name = "Alignment\nrate type", labels = c("Overall", "Unique"),
      begin = 0.2, end = 0.8, alpha = 0.5
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 15)
    )

  pdf(file.path(ratesDir, newFilename),
    width = 6, height = 4
  )
  print(boxplot)
  dev.off()
}
MakeBoxplot(ratesDf, "hisatBoxplot.pdf")
MakeBoxplot(ratesDf.naiveDrg, "hisatBoxplot.naiveDrg.pdf")
