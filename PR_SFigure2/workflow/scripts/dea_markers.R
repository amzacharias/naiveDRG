#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: scRNA of the naive human DRG from Renthal lab
# Author: Amanda Zacharias
# Date: 2025-07-22
# Email: 16amz1@queensu.ca
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.5.0

# Options -----------------------------------------------
setwd("/global/project/hpcg1553/Amanda_Zacharias/human_naive_drg")

# Packages -----------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(presto)

# Source -----------------------------------------------
source("workflow/scripts/helpers/theme_bw_bold.R")
source("workflow/scripts/helpers/save_plot.R")

# Pathways -----------------------------------------------
results_dir <- file.path("results")
## Input ===========
rdata_dir <- file.path(results_dir, "rdata")
neuron_data_human_path <- file.path(rdata_dir, "neuron_data_human_RNA.rds")

## Output ===========
plots_dir <- file.path(results_dir, "plots")
tables_dir <- file.path(results_dir, "tables")

# Load data -----------------------------------------------
neuron_data_human <- readRDS(neuron_data_human_path)

# Between cell-types DEA -----------------------------------------------
cat("\nRunning DEA\n")
## Run DEA ===========
human_markers <- FindAllMarkers(
  object = neuron_data_human,
  layer = "data",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.pct = 0.01,
  random.seed = 1
)

human_markers_clcn2 <- FindAllMarkers(
  object = neuron_data_human,
  features = "Clcn2",
  layer = "data",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.pct = 0.01,
  random.seed = 1
)
rownames(human_markers_clcn2) <- NULL
# Warning message:
# In FindMarkers.default(object = data.use, cells.1 = cells.1, cells.2 = cells.2,  :
#   No features pass logfc.threshold threshold; returning empty data.frame

human_markers_scn10a <- FindAllMarkers(
  object = neuron_data_human,
  layer = "data",
  features = "Scn10a",
  logfc.threshold = 0,
  test.use = "wilcox",
  min.pct = 0.01,
  random.seed = 1
)
rownames(human_markers_scn10a) <- NULL

## Multiple testing correction ===========
human_markers_clcn2_adj <- human_markers_clcn2 |>
  rename(p_val_adj_all_genes = p_val_adj) |>
  mutate(p_val_adj = p.adjust(p_val, method = "bonferroni", n = nrow(human_markers_clcn2)))

human_markers_scn10a_adj <- human_markers_scn10a |>
  rename(p_val_adj_all_genes = p_val_adj) |>
  mutate(p_val_adj = p.adjust(p_val, method = "bonferroni", n = nrow(human_markers_scn10a)))

## Subset results ===========
sig_human_markers <- human_markers |> subset(p_val_adj < 0.05)
human_markers_clcn2_adj |> subset(p_val_adj < 0.05)
human_markers_scn10a_adj |> subset(p_val_adj < 0.05)

## Save ===========
cat("\nSaving csv files\n")
write.csv(human_markers, file.path(tables_dir, "human_markers_all_res.csv"))
write.csv(sig_human_markers, file.path(tables_dir, "human_markers_all_sig_res.csv"))
write.csv(human_markers_clcn2_adj, file.path(tables_dir, "human_markers_clcn2_res.csv"))
write.csv(human_markers_scn10a_adj, file.path(tables_dir, "human_markers_scn10a_res.csv"))

# Save image -----------------------------------------------
cat("\nSaving image\n")
save.image(file = file.path(rdata_dir, "dea_markers.rData"))

# Session info -----------------------------------------------
sessionInfo()
