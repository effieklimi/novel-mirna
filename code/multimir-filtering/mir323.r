

library("ComplexHeatmap")
library("tidyverse")
library(plyr)


options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

miRnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply("as_tibble") %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) 
names(multimir) <- miRnames

mir323 <- multimir[[2]]


countFpkm <- read.csv(
  "results/tables/vsmcFpkm.csv",
  header = TRUE
) %>% as_tibble()

genes <- join(mir323[,1], countFpkm[, c(2,17:23)], by = "name") %>% na.omit

matrix <- genes[,2:7] 
rownames(matrix) <- make.names(genes$name, unique = TRUE)
matrix <- as.matrix(matrix) %>% +1 %>% log2()

mean <- rowMeans(matrix)
sd <- apply(matrix, 1, sd)
matrix <- (matrix-mean)/sd
colnames(matrix) <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p")




myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(500/2) + 1), 
              seq(max(matrix)/500, max(matrix), length.out = floor(500/2)))

pdf("results/figures/multimir-heatmap-mir323-all.pdf", width = 4, height = 10)
ComplexHeatmap::pheatmap(as.matrix(matrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    breaks=myBreaks
)
dev.off()