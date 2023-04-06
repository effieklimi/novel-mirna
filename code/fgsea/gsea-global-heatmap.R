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

miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

bp <- readRDS("results/rds/pathways/pathways-all-bioprocess.rds")
ke <- readRDS("results/rds/pathways/pathways-all-kegg.rds")
re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")

fgseaResultsSig <- read.csv(
    "results/tables/fgsea-vsmc-sigNES.csv",
    header = TRUE,
    row.names = 1
  )

heatmapGroupA <- fgseaResultsSig[complete.cases(fgseaResultsSig),]
heatmapGroupB <- fgseaResultsSig[rowSums(fgseaResultsSig, na.rm = TRUE) < 0 & rowSums(is.na(fgseaResultsSig)), ]
heatmapGroupC <- fgseaResultsSig[rowSums(fgseaResultsSig, na.rm = TRUE) > 0 & rowSums(is.na(fgseaResultsSig)), ]

fgseaHeatmap <- do.call("rbind", list(heatmapGroupA, heatmapGroupB, heatmapGroupC))

fgseaRankedGroupA <- heatmapGroupA[order(rowSums(heatmapGroupA), decreasing = FALSE),]
fgseaRankedGroupB <- heatmapGroupB[order(rowSums(heatmapGroupB, na.rm = TRUE), decreasing = FALSE),]
fgseaRankedGroupC <- heatmapGroupC[order(rowSums(heatmapGroupC, na.rm = TRUE), decreasing = FALSE),]

matchesGroupA <- which(rownames(fgseaHeatmap) %in% rownames(fgseaRankedGroupA[1:24,]))
matchesGroupB <- which(rownames(fgseaHeatmap) %in% rownames(fgseaRankedGroupB[1:19,]))
matchesGroupC <- which(rownames(fgseaHeatmap) %in% rownames(fgseaRankedGroupC[1:13,]))

matchesAll <- c(matchesGroupA, matchesGroupB, matchesGroupC)

colnames(fgseaHeatmap) <- miRNAnames

fgseaHeatmap <- fgseaHeatmap %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "GO BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))



assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
anot <- data.frame(row.names = rownames(fgseaHeatmap), Database = fgseaHeatmap[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "GO BP" = "grey90"))
ha <- rowAnnotation(foo = anno_mark(at = matchesAll, labels = rownames(fgseaHeatmap)[matchesAll]))

fgseaHeatmap[is.na(fgseaHeatmap)] <- 0

pdf(file = "results/figures/fgsea-vsmc-heatmap-global.pdf", width = 10, height = 15)
ComplexHeatmap::pheatmap(as.matrix(fgseaHeatmap[,c(2, 3, 5, 7, 1, 4, 6)]),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    right_annotation = ha,
    #na_col = "grey",
    annotation_row = anot,
    row_names_max_width = unit(18, "cm"),
    annotation_colors = ann_colors,
    row_split = rep(c(
      "A",
      "B",
      "C"), c(47, 38, 25)),
    colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaHeatmap[,c(2,3,5,7,1,4,6)] == 0, "·", ""), nrow(fgseaHeatmap[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()





































fgseaRanked <- as.matrix(fgseaResultsSig[order(rowSums(fgseaResultsSig), decreasing = FALSE),])
matches <- which(rownames(fgseaResultsSig) %in% rownames(fgseaRanked[1:10,]))












fgseaResultsSig <- fgseaResultsSig[order(rowSums(is.na(fgseaResultsSig))), ]


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

