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

# Whether Scn10a is a marker for the cell type ------------------------------------
# Scn10a is a marker for Nav1.8+ neurons
cells_scn10a_marker <- c(
  "Pvalb", "Mrgprd", "Ntrk3high+Ntrk2",
  "Calca+Sstr2", "Calca+Smr2", "Trpm8",
  "Ntrk3low+Ntrk2", "Atf3"
)
non_cells_scn10a_marker <- unique(
  neuron_data_human$Atlas_annotation
)[!unique(neuron_data_human$Atlas_annotation) %in% cells_scn10a_marker]

neuron_data_human$Nav1.8_marker <- ifelse(
  neuron_data_human$Atlas_annotation %in% cells_scn10a_marker,
  "Nav1.8_marker",
  "Not_Nav1.8_marker"
)
Idents(neuron_data_human) <- "Nav1.8_marker"
table(Idents(neuron_data_human))

# Between Nav1.8 DEA -----------------------------------------------
## DEA ===========
cat("\nRunning DEA\n")
dea_res <- FindMarkers(
  object = neuron_data_human,
  layer = "data",
  ident.1 = "Nav1.8_marker",
  ident.2 = "Not_Nav1.8_marker",
  test.use = "wilcox",
  features = "Clcn2"
)

## Correct only for gene tested, not all genes in the dataset ===========
dea_res_adj <- dea_res |>
  rename(p_val_adj_all_genes = p_val_adj) |>
  mutate(p_val_adj = p.adjust(p_val, method = "bonferroni", n = nrow(dea_res)))

dea_res_sig <- dea_res_adj |> subset(p_val_adj < 0.05)

## Save ===========
cat("\nSaving csv files\n")
write.csv(dea_res_adj, file.path(tables_dir, "dea_nav1.8_human_res.csv"))
write.csv(dea_res_sig, file.path(tables_dir, "dea_nav1.8_human_sig_res.csv"))

# Plotting of dea results -----------------------------------------------
cat("\nPlotting\n")
plot_data <- FetchData(
  object = neuron_data_human,
  vars = c("Clcn2", "Atlas_annotation", "Nav1.8_marker")
)
violin_dea_plot <- plot_data |>
  arrange(Nav1.8_marker, Atlas_annotation) |>
  mutate(Atlas_annotation = factor(Atlas_annotation, levels = unique(Atlas_annotation))) |>
  ggplot(aes(x = Atlas_annotation, y = Clcn2, fill = Nav1.8_marker)) +
  geom_violin(alpha = 0.5, scale = "width", adjust = 0.2) +
  geom_jitter(
    aes(color = Nav1.8_marker),
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
  ) +
  scale_y_continuous(
    name = expression(italic("CLCN2")),
    limits = c(0, NA), expand = c(0.01, 0.01)
  ) +
  scale_fill_manual(
    values = c("Nav1.8_marker" = "#00BFC4", "Not_Nav1.8_marker" = "#F8766D"),
    labels = c("Nav1.8_marker" = "yes", "Not_Nav1.8_marker" = "no")
  ) +
  scale_color_manual(
    values = c("Nav1.8_marker" = "#00BFC4", "Not_Nav1.8_marker" = "#F8766D"),
    labels = c("Nav1.8_marker" = "yes", "Not_Nav1.8_marker" = "no")
  ) +
  theme_bw_bold() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
save_plot(
  plot_obj = violin_dea_plot,
  new_filename = "violin_dea_nav1.8_clcn2.pdf",
  new_path = plots_dir,
  w = 360,
  h = 180
)

# Save image -----------------------------------------------
cat("\nSaving image\n")
save.image(file = file.path(rdata_dir, "dea_nav1.8.rData"))

# Session info -----------------------------------------------
sessionInfo()
