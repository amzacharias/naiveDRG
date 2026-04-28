#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Custom Theme adapted from theme_bw()
# Author: Amanda Zacharias
# Date: 2025-06-06
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(ggplot2) # 3.5.1

# Source -----------------------------------------------

# Pathways -----------------------------------------------
## Input ===========

## Output ===========

# Load data -----------------------------------------------
theme_bw_bold <- function() {
  theme_bw() %+replace%
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = rel(1.2), hjust = 0.5),
      axis.title = element_text(size = rel(1)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.25),
      axis.ticks = element_line(linewidth = 1),
      axis.text = element_text(colour = "black"),
      plot.caption = element_text(size = rel(1))
    )
}
