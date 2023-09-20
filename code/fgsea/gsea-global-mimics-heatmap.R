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
names(bp) %>% lapply(str_to_sentence) %>% lapply(gsub, pattern = "Dna", replacement = "DNA", fixed = TRUE) -> names(bp)

fgseaResultsSig <- read.csv(
    "results/tables/fgsea-vsmc-sigNES-mimics.csv",
    header = TRUE,
    row.names = 1
  )


# Pathways affected by all miRNAs: # ~ NES no NA:
fgseaHeatmapNesNeg <- fgseaResultsSig[complete.cases(fgseaResultsSig), ]  %>% write.csv("results/tables/heatmapGroupA.csv")
# Pathways affected by at least 1 miRNA but not all: # ~ NES with NA:
fgseaHeatmapNesPos <- fgseaResultsSig[!complete.cases(fgseaResultsSig), ]

fgseaRankedGroupA <- fgseaHeatmapNesNeg[order(rowSums(fgseaHeatmapNesNeg), decreasing = FALSE), ] %>% write.csv("results/tables/heatmapGroupA.csv")
fgseaRankedGroupB <- fgseaHeatmapNesPos[order(rowSums(fgseaHeatmapNesPos, na.rm = FALSE), decreasing = FALSE), ] %>% write.csv("results/tables/heatmapGroupB.csv")
fgseaRankedGroupC <- heatmapGroupC[order(rowSums(heatmapGroupC, na.rm = TRUE), decreasing = TRUE), ] %>% write.csv("results/tables/heatmapGroupC.csv")

matchesNesNeg <- which(rownames(fgseaHeatmapNesNeg) %in% rownames(fgseaRankedGroupA[1:(nrow(fgseaRankedGroupA) / 3), ]))
matchesNesPos <- which(rownames(fgseaHeatmapNesNeg) %in% rownames(fgseaRankedGroupB[1:(nrow(fgseaRankedGroupB) / 3), ]))
#matchesGroupC <- which(rownames(fgseaHeatmapNesPos) %in% rownames(fgseaRankedGroupC[1:(nrow(fgseaRankedGroupC) / 2), ]))

colnames(fgseaHeatmapNesNeg) <- miRNAnames
colnames(fgseaHeatmapNesPos) <- miRNAnames

fgseaHeatmapNesNeg <- fgseaHeatmapNesNeg %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "GO BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))

fgseaHeatmapNesPos <- fgseaHeatmapNesPos %>%
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
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "GO BP" = "grey90"))

anotNesNeg <- data.frame(row.names = rownames(fgseaHeatmapNesNeg), Database = fgseaHeatmapNesNeg[, 8])
#haNesNeg <- rowAnnotation(foo = anno_mark(at = matchesNesNeg, labels = rownames(fgseaHeatmapNesNeg)[matchesNesNeg]))

haNesNeg <- rowAnnotation(
  foo1 = anno_mark(
    at = matchesNesNeg,
    labels = rownames(fgseaHeatmapNesNeg)[matchesNesNeg],
    labels_gp = gpar(fontsize = 13),
    link_gp = gpar(lwd = 0.7)
    )
  )

fgseaHeatmapNesNeg[is.na(fgseaHeatmapNesNeg)] <- 0
myBreaks <- c(seq(min(as.matrix(fgseaHeatmapNesNeg[, c(2, 3, 5, 7, 1, 4, 6)])), 0, length.out=ceiling(500/2) + 1), 
              seq(max(as.matrix(fgseaHeatmapNesNeg[, c(2, 3, 5, 7, 1, 4, 6)]))/500, max(as.matrix(fgseaHeatmapNesNeg[, c(2, 3, 5, 7, 1, 4, 6)])), length.out = floor(500/2)))


pdf(file = "results/figures/fgsea-vsmc-heatmap-global.pdf", width = 12, height = 15)
ComplexHeatmap::pheatmap(as.matrix(fgseaHeatmapNesNeg[, c(2, 3, 5, 7, 1, 4, 6)]),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 5,
    fontsize_col = 9,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    right_annotation = haNesNeg,
    #na_col = "grey",
    annotation_row = anotNesNeg,
    row_names_max_width = unit(18, "cm"),
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    border_gp = gpar(col = "black", lwd = 2),
    breaks = myBreaks,
    display_numbers = matrix(ifelse(fgseaHeatmapNesNeg[,c(2,3,5,7,1,4,6)] == 0, "·", ""), nrow(fgseaHeatmapNesNeg[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()



anotNesPos <- data.frame(row.names = rownames(fgseaHeatmapNesPos), Database = fgseaHeatmapNesPos[, 8])

haNegPos <- rowAnnotation(
  foo1 = anno_mark(
    at = matchesNesPos,
    labels = rownames(fgseaHeatmapNesPos)[matchesNesPos],
    labels_gp = gpar(fontsize = 16),
    link_gp = gpar(lwd = 0.7)
    )
  )

fgseaHeatmapNesPos[is.na(fgseaHeatmapNesPos)] <- 0
myBreaks <- c(seq(min(as.matrix(fgseaHeatmapNesPos[, c(2, 3, 5, 7, 1, 4, 6)])), 0, length.out=ceiling(500/2) + 1), 
              seq(max(as.matrix(fgseaHeatmapNesPos[, c(2, 3, 5, 7, 1, 4, 6)]))/500, max(as.matrix(fgseaHeatmapNesPos[, c(2, 3, 5, 7, 1, 4, 6)])), length.out = floor(500/2)))



pdf(file = "results/figures/fgsea-vsmc-heatmap-global.pdf", width = 10, height = 13)
ComplexHeatmap::pheatmap(as.matrix(fgseaHeatmapNesPos[, c(2, 3, 5, 7, 1, 4, 6)]),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 5,
    fontsize_col = 15,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    right_annotation = haNegPos,
    #na_col = "grey",
    annotation_row = anotNesPos,
    row_names_max_width = unit(18, "cm"),
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    breaks = myBreaks,
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaHeatmapNesPos[, c(2,3,5,7,1,4,6)] == 0, "·", ""), nrow(fgseaHeatmapNesPos[,c(2,3,5,7,1,4,6)]))
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

