

library("ComplexHeatmap")
library("tidyverse")
library(plyr)

miRnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply("as_tibble") %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) 
names(multimir) <- miRnames
multimir <- multimir[c(2, 3, 5, 7, 1, 4, 6)]

targets <- do.call(rbind, multimir)

countFpkm <- read.csv(
  "results/tables/vsmcFpkm.csv",
  header = TRUE
) %>% as_tibble()

genes <- join(targets[,1], countFpkm[, c(2,17:40)], by = "name") %>% na.omit
genesDuplicates <- 
    split(genes, duplicated(genes$name) | duplicated(genes$name, fromLast = TRUE)) %>%
    lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
    lapply(data.frame) %>%
    lapply(column_to_rownames, var = "name")

matrixUnique <- genesDuplicates$`FALSE` %>% as.matrix() %>% +1 %>% log2()
mean <- rowMeans(matrixUnique)
sd <- apply(matrixUnique, 1, sd)
matrixUnique <- (matrixUnique-mean)/sd
colnames(matrixUnique) <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")

matrixDuplicate <- genesDuplicates$`TRUE` %>% as.matrix() %>% +1 %>% log2()
mean <- rowMeans(matrixDuplicate)
sd <- apply(matrixDuplicate, 1, sd)
matrixDuplicate <- (matrixDuplicate-mean)/sd
colnames(matrixDuplicate) <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")

matrix <- genes[,2:25] 
rownames(matrix) <- make.names(genes$name, unique = TRUE)
matrix <- as.matrix(matrix) %>% +1 %>% log2()

mean <- rowMeans(matrix)
sd <- apply(matrix, 1, sd)
matrix <- (matrix-mean)/sd
colnames(matrix) <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")








myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(500/2) + 1), 
              seq(max(matrix)/500, max(matrix), length.out = floor(500/2)))

pdf("results/figures/multimir-heatmap-all.pdf", width = 7, height = 10)
ComplexHeatmap::pheatmap(as.matrix(matrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    row_split = rep(c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7"), c(492, 1066, 570, 960, 1271, 273, 246)),
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    breaks=myBreaks
)
dev.off()








myBreaks <- c(seq(min(matrixUnique), 0, length.out=ceiling(500/2) + 1), 
              seq(max(matrixUnique)/500, max(matrixUnique), length.out = floor(500/2)))

pdf("results/figures/multimir-heatmap-unique.pdf", width = 6, height = 10)
ComplexHeatmap::pheatmap(as.matrix(matrixUnique),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    row_split = rep(c(
    "miR-323a-3p",
    "miR-449b-5p",
    "miR-491-3p",
    "miR-892b",
    "miR-1827",
    "miR-4774-3p",
    "miR-5681b"), c(482, 1056, 560, 950, 1261, 263, 228)),
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    breaks=myBreaks
)
dev.off()




myBreaks <- c(seq(min(matrixDuplicate), 0, length.out=ceiling(500/2) + 1), 
              seq(max(matrixDuplicate)/500, max(matrixDuplicate), length.out = floor(500/2)))

pdf("results/figures/multimir-heatmap-duplicate.pdf", width = 6, height = 7)
ComplexHeatmap::pheatmap(as.matrix(matrixDuplicate),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    row_split = rep(c(
    "miR-323a-3p",
    "miR-449b-5p",
    "miR-491-3p",
    "miR-892b",
    "miR-1827",
    "miR-4774-3p",
    "miR-5681b"), c(482, 1056, 560, 950, 1261, 263, 228)),
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500),
    breaks=myBreaks
)
dev.off()

