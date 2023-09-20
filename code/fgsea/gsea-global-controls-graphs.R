#!/usr/bin/env Rscript --vanilla
library('ComplexHeatmap')
library('tidyverse')

##############################

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

########## pathways ##########
draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)
}

controlNames <- c("IL-1a/PDGF-BB", "Mock", "miR-CTRL")
controlNames <- c("pathway", "ip", "Mock", "miRCTRL")

bp <- readRDS("results/rds/pathways/pathways-all-bioprocess.rds")
ke <- readRDS("results/rds/pathways/pathways-all-kegg.rds")
re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")
names(bp) %>% lapply(str_to_sentence) %>% lapply(gsub, pattern = "Dna", replacement = "DNA", fixed = TRUE) -> names(bp)

fgseaResults <- read.csv(
    "results/tables/fgsea-vsmc-sigNES-controls.csv",
    header = TRUE,
    row.names = 1
  ) %>% rownames_to_column(var = "pathway") %>% as_tibble() %>% na.omit() 

colnames(fgseaResults) <- controlNames

# Pathways affected by at least 1 miRNA but not all, < 0 NES:
heatmapGroupNeg <- fgseaResults[fgseaResults$miRCTRL < 0, ]
heatmapGroupPos <- fgseaResults[fgseaResults$miRCTRL > 0, ]
heatmapGroupNeg <- heatmapGroupNeg[order(heatmapGroupNeg$miRCTRL), ]
heatmapGroupPos <- heatmapGroupPos[order(heatmapGroupPos$miRCTRL, decreasing = TRUE), ]

pdf("results/figures/gsea-ipCTRL-negativeNes.pdf", width = 8, height = 4)
heatmapGroupNeg %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>% 
  dplyr::arrange(-dplyr::row_number()) %>%
  ggplot(aes(y = miRCTRL, x = factor(pathway, levels = pathway))) +
      geom_bar(stat = "identity", fill = "#3c73d1", alpha = 1, width = .7) +
      theme_minimal() +
      xlab("") +
      ylab("Normalised Enrichment Score") +
      scale_y_reverse() +
      coord_flip()
dev.off()

pdf("results/figures/gsea-ipCTRL-positiveNes.pdf", width = 8, height = 9)
heatmapGroupPos[c(1:60), ] %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>% 
  dplyr::arrange(-dplyr::row_number()) %>%
  ggplot(aes(x = miRCTRL, y = factor(pathway, levels = pathway))) +
      geom_bar(stat = "identity", fill = "#c13b89", alpha = .6, width = .7) +
      theme_minimal() +
      xlab("Normalised Enrichment Score") +
      ylab("")
dev.off()

heatmapGroupPos[c(1:60), ] %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>% 
  write.csv("results/tables/fgsea-ipCTRL-positiveNes.csv")