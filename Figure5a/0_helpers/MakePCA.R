#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make PCA with DESeq2 helper functions
# Author: Amanda Zacharias
# Date: 2023-11-30
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(DESeq2) # 1.26.0


# GetPcaData -----------------------------------------
GetPcaData <- function(vstDat, removeVec, toPlotVars = c("sampleGrps")) {
  #' This function performs principal component analysis and stores results
  #' in a way that is compatible with ggplot2 for visualization.
  #' @param vstData Large DESeqTransform object with VS transformed counts
  #' @param removeVec Sample number of any sample that is an outlier.
  #' Assumes the 'sampleNum' column from metadata has this information!!
  #' @param toPlotVars Character vector with the variables to be included 
  #' in output for plotting.
  #' @return Dataframe fit for ggplot2
  lbls <- rep(FALSE, ncol(vstDat))
  lbls[match(removeVec, vstDat$sampleNum)] <- TRUE
  pcaData <- plotPCA(vstDat, intgroup = toPlotVars, returnData = TRUE) %>%
    mutate("labels" = lbls)
  return(pcaData)
}

# MakePCA -----------------------------------------
MakePCA <- function(pcaDat, toPlotVars = c("sampleGrps", "ztTime", "labels"),
                    newFilename = "pca", newPath = ".") {
  #' Generates a plot for PCA analysis and save as a pdf
  percentVar <- round(100 * attr(pcaDat, "percentVar"))
  for (group in toPlotVars) {
    if (sum(pcaDat$labels) >= 1 & group != "labels") {
      pcaGplot <- pcaDat %>%
        ggplot(aes(PC1, PC2, color = labels, fill = .data[[group]])) +
        geom_point(shape = 21, size = 5, alpha = 0.6, stroke = 1) +
        scale_color_manual(
          name = "Is outlier", labels = c("no", "yes"),
          values = c("white", "black")
        ) +
        guides(color = guide_legend(order = 1, byrow = TRUE), fill = guide_legend(order = 2))
    } else {
      pcaGplot <- pcaDat %>%
        ggplot(aes(PC1, PC2, color = .data[[group]])) +
        geom_point(size = 5, alpha = 0.6)
    }
    pcaGplot <- pcaGplot +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(str_to_sentence(group)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(-1, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.key = element_rect(fill = "grey90"),
        text = element_text(size = 15)
      )
    if (length(groups) > 6) {
      pcaGplot <- theme(legend.position = "none")
    }
    ggsave(
      plot = pcaGplot, path = newPath,
      filename = paste(group, newFilename, sep = "."),
      width = 185, height = 185, units = "mm", dpi = 300
    )
  }
  cat("\nSaved to:", newPath)
}


