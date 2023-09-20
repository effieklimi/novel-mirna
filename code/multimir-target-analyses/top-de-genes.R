library("ComplexHeatmap")
library("tidyverse")
library("org.Hs.eg.db")
library("enrichR")
library("clusterProfiler")
library(dendextend)
library(plyr)

#
options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
# naming lists
miRnames <- c("hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-491-3p", "hsa-miR-892b", "hsa-miR-1827", "hsa-miR-4774-3p", "hsa-miR-5681b")
colNames <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")
#

multimir <-
  readRDS("results/rds/vsmc-multimir.rds")[c(2, 3, 5, 7, 1, 4, 6)] %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply("[", c(1:20), )
names(multimir) <- miRnames
multimir <- multimir[c(2, 3, 5, 7, 1, 4, 6)]
targets <- do.call(rbind, multimir)

genes <- join(targets[, 1], countFpkm[, c(2, 17:40)], by = "name")
matrix <- genes[, 2:25]
rownames(matrix) <- make.names(genes$name, unique = TRUE)
matrix <- as.matrix(matrix) %>% +1 %>% log2()
mean <- rowMeans(matrix)
sd <- apply(matrix, 1, sd)
matrix <- (matrix - mean) / sd
colnames(matrix) <- colNames
#
pdf("results/figures/multimir-heatmap-topgenes.pdf", width = 7, height = 10)
ComplexHeatmap::pheatmap(as.matrix(matrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 4,
    fontsize_col = 7,
    row_split = rep(c(
      "1",
      "2",
      "3",
      "4",
      "5",
      "6",
      "7"), c(20, 20, 20, 20, 20, 20, 20), levels(miRnames)),
    row_title_gp = gpar(fontsize = 8),
    show_rownames = TRUE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500)
)
dev.off()