#!/usr/bin/env Rscript --vanilla
library('ComplexHeatmap')
library('tidyverse')

##############################

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

########## pathways ##########
read.geneset <- function(path_to_gset)  {
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]] != ""])
  names(out) = bp$geneset.names
  return(out)
}

bp <- readRDS("results/rds/pathways/pathways-all-bioprocess.rds")
ke <- readRDS("results/rds/pathways/pathways-all-kegg.rds")
re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")

draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)
}

fgseaResultsSig <- read.csv(
    "results/tables/fgsea-vsmc-sigNES-p01.csv",
    header = TRUE, 
    row.names = 1
  )

fgseaResultsSig[is.na(fgseaResultsSig)] <- 0

# Rank results based on decreasing number of NAs per row and format for heatmap:
fgseaRanked <- as.matrix(fgseaResultsSig[order(rowSums(fgseaResultsSig), decreasing = FALSE),])
matches <- which(rownames(fgseaResultsSig) %in% rownames(fgseaRanked[1:40,]))
fgseaMatrix <- as.matrix(fgseaResultsSig)
#rownames(fgseaMatrix) <- rownames(fgseaResultsSig)
fgseaMatrix <- fgseaMatrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))


assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
anot <- data.frame(row.names = rownames(fgseaResultsSig), Database = fgseaMatrix[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "Gene Ontology BP" = "grey90"))
ha <- rowAnnotation(foo = anno_mark(at = matches, labels = rownames(fgseaMatrix)[matches]))

pdf(file = "results/figures/fgsea-vsmc-heatmap-p01.pdf", width = 10, height = 13)
ComplexHeatmap::pheatmap(as.matrix(fgseaMatrix[,c(2, 3, 5, 7, 1, 4, 6)]),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 7,
    fontsize_col = 7,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "grey",
    annotation_row = anot,
    row_names_max_width = unit(12, "cm"),
    annotation_colors = ann_colors,
    row_split = rep(c("A", "B", "C"), c(26, 26, 68)),
    main = "Gene Set Enrichment Analysis",
    colorRampPalette(c("#005dc7", "white", "#d60047"))(300),
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaMatrix[,c(2,3,5,7,1,4,6)] == 0, "·", ""), nrow(fgseaMatrix[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()