#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make volcano plot helper
# Author: Amanda Zacharias
# Date: 2023-04-09
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
options(ggrepel.max.overlaps = Inf, box.padding = 0)

# Packages -----------------------------------------
library(ggplot2) # 3.4.1
library(ggrepel) # 0.9.3

# Pathways -----------------------------------------



# Load data -----------------------------------------



# Function -----------------------------------------
MakeVolcano <- function(df, showName = TRUE,
                        newTitle = "Volcano", newPath = ".", fileEnd = "pdf",
                        idColumn = "isoform_id",
                        newWidth = 185, newHeight = 185,
                        showTitle = TRUE,
                        toMoveLegend = FALSE
                        ) {
  #' Make a volcano plot
  #'
  # Copy results table to new variable that we can alter
  volcanoDf <- df
  # Mark significant rows
  volcanoDf$sig <- volcanoDf$padj < 0.05
  # Order dataframe and get top 20 sig transcripts'
  if (showName == TRUE) {
    num2label <- 15
  } else if (showName == FALSE) {
    # transcript ids take up more space, so label less
    num2label <- 10
  }
  resOrdered <- volcanoDf %>% arrange(pvalue)
  resOrdered$labels <- resOrdered[, idColumn] %in%
    resOrdered[1:num2label, idColumn]

  # Plotting
  bcThreshold <- 0.05 / nrow(resOrdered)
  cat("\nBonferroni threshold is:", bcThreshold, "\n")

  volcano <- resOrdered %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
    annotate("rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = Inf,   fill = "#ffffff00", alpha = 0.2) +  # #CDB293
    annotate("rect", xmin = 0, xmax = Inf,   ymin = -Inf, ymax = Inf,   fill = "#939393", alpha = 0.2) +  # "#386FA3
    geom_point(aes(colour = sig), size = 2.5, show.legend = FALSE) + # size = 1.75
    ggtitle(newTitle) +
    xlab(expression("Log"[2] * "(fold change)")) +
    ylab(expression("-Log"[10] * " (p-value)")) +
    geom_hline(aes(yintercept = -log10(bcThreshold), linetype = "Adjusted"), linewidth = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), linetype = "Nominal"), linewidth = 0.8) +
    scale_linetype_manual(name = "Significance", values = c(2, 4)) +
    scale_color_manual(
      values = c("#000000", "#FF0000"),
      labels = c("Not significant", "Significant")
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25)),
      axis.text = element_text(colour = "#000000")
    )
  # Labels
  if (showName == TRUE) {
    volcano <- volcano + geom_label_repel(
      # size = 5,
      size = 3, nudge_y = 2, min.segment.length = 0.25,
      force = 10, show.legend = FALSE,
      aes(label = ifelse(labels == TRUE, # label top t's.
        as.character(gene_name), ""
      ))
    )
  } else if (showName == FALSE) {
    volcano <- volcano + geom_label_repel(
      # size = 5,
      size = 3, nudge_y = 2, min.segment.length = 0.25,
      force = 10, show.legend = FALSE,
      aes(label = ifelse(labels == TRUE, # label top t's.
        as.character(.data[[idColumn]]), ""
      ))
    )
  }
  if (showTitle == FALSE) {
    volcano <- volcano + theme(plot.title = element_blank())
  }
  if (toMoveLegend == TRUE) {
    volcano <- volcano + theme(legend.position = c(0.22, 0.85))
  }
  # Save
  ggsave(
    plot = volcano, path = newPath,
    filename = paste("showname", showName, "volcano", fileEnd, sep = "."),
    width = newWidth, height = newHeight, units = "mm"
  )
}
