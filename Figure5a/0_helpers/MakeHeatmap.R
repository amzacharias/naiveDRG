#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make heatmap gene counts across samples
# Author: Amanda Zacharias
# Date: 2023-07-13
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------

# MakeSampleHeatmap -----------------------------------------
MakeSampleHeatmap <-function(countMat, varOfInterest, coldat, newFilename, newPath) {
  #' Custom function to make a heatmap of samples' distance
  #'
  #' @param countMat gene count matrix
  #' @param varOfInterest variable of interest (string); matches in coldat colnames
  #' @param coldat metadata for samples
  #' @param newFilename filename for output .pdf file (string)
  #' @param newPath path to output folder (string)
  #' @return Creates a pdf file
  #' @example
  #'
  # Distance calculation
  distObj <- dist(t(countMat))
  distMat <- distObj %>% as.matrix()
  rownames(distMat) <- coldat[[varOfInterest]]
  colnames(distMat) <- NULL
  # Make plot
  hplot <- pheatmap(distMat,
                    color = cividis(n = 100),
                    show_rownames = TRUE,
                    clustering_distance_rows = distObj,
                    clustering_distance_cols = distObj
  )
  # Save
  ggsave(
    plot = hplot, filename = newFilename, path = newPath,
    width = 92.5, height = 92.5, units = "mm"
  )
  cat("\nPlot saved to:", newPath)
} 

# MakeGeneHeatmap -----------------------------------------
# Column annotation colours =================
# 0-9 = light, 12-24 = dark
# drkColors <- colorRampPalette(c("gray50", "gray15"))(2)
# lgtColors <- colorRampPalette(c("lightgoldenrod", "orange2"))(2)
# annoColors <- c(lgtColors, drkColors)
# names(annoColors) <-
#   as.factor(unique(coldata$ztTime))
# annoColors <- list(time = annoColors)


MakeGeneHeatmap <-
  function(countMat,
           varOfInterest, 
           coldat,
           title,
           newFilename,
           newPath,
           annoCol,
           cols = cividis(150),
           brks = NA,
           dim = list(width = 4, height = 7),
           showrownames = FALSE,
           clustRows = FALSE,
           newfontsize = 10) {
    #' Custom function to make an expression heatmap for genes/transcripts
    #' Assumes that order of samples in count matrix and coldat is the same!
    #'
    #' @param countMat Input count matrix (likely VST counts)
    #' @param varOfInterest Variable that you want to annotate by.
    #' @param coldat Metadata for samples (dataframe)
    #' @param title Title for plot (string)
    #' @param newFilename Filename prefix (string)
    #' @param newPath Path to output directory (string)
    #' @param annoCol What colours should the samples be annotated with (list)
    #' @param brks Breakslist for color range (allow same legend across plots)
    #' @param cols Vector colors
    #' @param dim Dimensions of output file
    #' @param showrownames  Whether to show rownames in plot (TRUE/FALSE)
    #' @param clustRows Whether to cluster shows by euclidian distance
    #' @param newfontsize Font size for text in plot
    #' @return Writes a pdf file.
    #' @example
    #' MakeHeatmap(df, "ztTime", coldata, "title", "genes", "./plots", annoColors)
    # Colors
    annoDf <- data.frame(
      time = as.factor(coldat[[varOfInterest]]),
      row.names = colnames(countMat)
    )
    # Plot
    hPlot <- pheatmap(
      mat = countMat, 
      color = cols,
      breaks = brks,
      border_color = NA,
      scale = "row",
      show_rownames = showrownames,
      show_colnames = FALSE,
      annotation_col = annoDf,
      annotation_colors = annoCol,
      annotation_names_col = FALSE,
      fontsize = newfontsize,
      fontsize_row = 10,
      fontsize_col = 10,
      main = title,
      angle_col = 45,
      cluster_rows = clustRows,
      cluster_cols = FALSE,
      annotation_legend = TRUE,
      height = dim[["height"]],
      width = dum[["width"]]
    )
    ggsave(
      plot = hPlot,
      filename = newFilename,
      path = newPath,
      width = dim[["width"]],
      height = dim[["height"]]
    )
    cat("\nPlot saved to:", newPath)
  }

