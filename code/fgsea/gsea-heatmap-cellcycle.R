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

fgseaResCellCycle <- read.csv(
  "results/tables/fgsea-vsmc-cellcycle-sigNES.csv",
  header = TRUE,
  row.names = 1
)

fgseaResCellCycle <- fgseaResCellCycle[order(rowSums(is.na(fgseaResCellCycle))), ]

fgseaResCellCycle[is.na(fgseaResCellCycle)] <- 0
fgseaRankedCellCycle <- as.matrix(fgseaResCellCycle[order(rowSums(fgseaResCellCycle), decreasing = FALSE),])
#matchesCellCycle <- which(rownames(fgseaResCellCycle) %in% rownames(fgseaRankedCellCycle[1:10,]))
fgseaMatrixCellCycle <- as.matrix(fgseaResCellCycle)
#rownames(fgseaMatrix) <- rownames(fgseaResCellCycle)
fgseaMatrixCellCycle <- fgseaMatrixCellCycle %>%
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
anot <- data.frame(row.names = rownames(fgseaResCellCycle), Database = fgseaMatrixCellCycle[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "Gene Ontology BP" = "grey90"))
#ha <- rowAnnotation(foo = anno_mark(at = matchesCellCycle, labels = rownames(fgseaMatrixCellCycle)[matches]))
pdf(file = "results/figures/fgsea-vsmc-heatmap-cellcycle.pdf", width = 11, height = 8)
ComplexHeatmap::pheatmap(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)],
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 12,
    fontsize_col = 12,
    show_rownames = TRUE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    annotation_row = anot,
    row_names_max_width = unit(18, "cm"),
    annotation_colors = ann_colors,
    colorRampPalette(c("#002f80", "white"))(300),
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)] == 0, "·", ""), nrow(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()