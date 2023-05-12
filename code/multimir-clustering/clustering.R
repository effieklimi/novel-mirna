library("ComplexHeatmap")
library("tidyverse")
library("org.Hs.eg.db")
library("enrichR")
library("clusterProfiler")
library(dendextend)
library(plyr)


options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")


# naming lists
miRnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
colNames <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")
#

# loading expression data and multimir data
multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .5)) # keeping the top 50% of genes per miRNA based on LFC
names(multimir) <- miRnames
multimir <- multimir[c(2, 3, 5, 7, 1, 4, 6)]
targets <- do.call(rbind, multimir)

countFpkm <- read.csv(
  "results/tables/vsmcFpkm.csv",
  header = TRUE
) %>% as_tibble()

genes <- join(targets[, 1], countFpkm[, c(2, 17:40)], by = "name") %>% na.omit %>% distinct(name, .keep_all = TRUE)
matrix <- genes[, 2:25]
rownames(matrix) <- genes$name
matrix <- as.matrix(matrix) %>% +1 %>% log2()
mean <- rowMeans(matrix)
sd <- apply(matrix, 1, sd)
matrix <- (matrix - mean) / sd
colnames(matrix) <- colNames
#


# clustering 
distances <- hclust(dist(matrix))
clusters <- cut_lower_fun(as.dendrogram(distances), h = 9.05)
#9.14 gives 6(114-467) mean = 377.5
#9.13 gives 6(114-467) mean = 377.5
#9.12 gives 7(114-467) median = 280
#9.11 gives 7(114-467) median = 280
#9.07  gives 8(114-467) median = 214
#9.05  gives 8(114-467) median = 214 <-
#9.04  gives 9(66-467) median = 188

length(clusters)
min(unlist(lapply(clusters, length)))
max(unlist(lapply(clusters, length)))
median(unlist(lapply(clusters, length)))
#


# making DF of gene numbers per cluster for heatmap breaks
clusterGeneNo <- unlist(lapply(clusters, length))
heatmapBreaks <- c()
heatmapBreaks[1] <- clusterGeneNo[1]
for (i in 2:length(clusterGeneNo)) {
    heatmapBreaks[i] <- heatmapBreaks[i - 1] + clusterGeneNo[i]
}
clustersNo <- data.frame(
    name = seq_along(clusters),
    geneNo = clusterGeneNo,
    heatmapBreaks = heatmapBreaks,
    stringsAsFactors = FALSE
)
#


# making of matrix for heatmap containing clusters
clusterExpr <- list()
for (i in seq_len(length(clusters))) {
    clusterExpr[[i]] <- filter(countFpkm, name %in% clusters[[i]])
}
clusterMatrices <- do.call(rbind, clusterExpr)
clusterMatrices <- clusterMatrices[,c(2, 17:40)] %>% na.omit %>% distinct(name, .keep_all = TRUE)
clusterMatrix <- as.data.frame(clusterMatrices[,2:25])
rownames(clusterMatrix) <- clusterMatrices$name
clusterMatrix <- as.matrix(clusterMatrix) %>% +1 %>% log2()
mean <- rowMeans(clusterMatrix)
sd <- apply(clusterMatrix, 1, sd)
clusterMatrix <- (clusterMatrix - mean) / sd
colnames(clusterMatrix) <- colNames
#



# heatmap
pdf("results/figures/multimir-heatmap-all-new.pdf", width = 7, height = 10)
ComplexHeatmap::pheatmap(as.matrix(clusterMatrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    row_split = rep(clustersNo$name, clustersNo$geneNo),
    legend = TRUE,
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500)
)
dev.off()
#


geneOntology <-
  clusters %>%
  map(enrichGO,
      universe      = countFpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      minGSSize     = 10)

geneOntologyTable <- lapply(geneOntology, function(x) x@result)

write.table(geneOntologyTable[[1]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster1.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[2]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster2.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[3]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster3.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[4]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster4.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[5]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster5.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[6]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster6.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[7]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster7.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneOntologyTable[[8]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/clusterData/cluster8.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))


"/Users/effieklimi/Documents/novel-mirna/GO-Figure/clusterData"














addOne <- function(x) {
    x + 1
}





myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(500/2) + 1), 
              seq(max(matrix)/500, max(matrix), length.out = floor(500/2)))

pdf("results/figures/multimir-heatmap-all-new.pdf", width = 7, height = 10)
ComplexHeatmap::pheatmap(as.matrix(matrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = distances,
    fontsize_row = 10,
    fontsize_col = 12,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    colorRampPalette(c("#002f80", "white", "#87003f"))(500)
)
dev.off()


clusterMatrices <- 
    clusterExpr %>% 
    lapply( "[", , c(2, 17:40)) %>%
    lapply(na.omit) %>%
    lapply(distinct, name, .keep_all = TRUE) %>%
    lapply(as.data.frame) %>%
    lapply(column_to_rownames, "name") %>%
    lapply(as.matrix) %>%
    lapply(addOne) %>%
    lapply(log2)


means <- lapply(clusterMatrices, rowMeans)
sd <- lapply(clusterMatrices, function(y) apply(y, 1, sd))
matrix <- (matrix-means)/sd

lapply(clusterMatrices, function(y) apply(y, 1, sd))




#9 gives 12(91-567) mean = 321
#9.13 gives 8(91-785) mean = 321
#9.135 gives 7(105-785)
#9.14 gives 7(105-785)
#9.15 gives 7(105-785)

# For the pathway analysis we took the uppermost quartile wrt logfoldchange
geneOnyology <-
  clusters %>%
  map(enrichGO,
      universe      = countFpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff = 1,
      minGSSize = 10)

bpEnrichTable <- lapply(geneOnyology, function(x) x@result)

bpEnrichTable[[1]][c(2,6),1:10]
bpEnrichTable[[2]][c(2,6),1:10]
bpEnrichTable[[3]][c(2,6),1:10]
bpEnrichTable[[4]][c(2,6),1:10]
bpEnrichTable[[5]][c(2,6),1:10]
bpEnrichTable[[6]][c(2,6),1:10]
bpEnrichTable[[7]][c(2,6),1:10]
bpEnrichTable[[8]][c(2,6),1:10]
bpEnrichTable[[9]][c(2,6),1:10]
bpEnrichTable[[10]][c(2,6),1:10]
bpEnrichTable[[11]][c(2,6),1:10]
bpEnrichTable[[12]][c(2,6),1:10]
bpEnrichTable[[13]][c(2,6),1:10]




multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)

genes <-
  multimir %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)
names(genes) <- miRNAnames

re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")
names(re) <- tolower(names(re))


fgseaRes <- 
  map(genes, fgsea, pathways = re, minSize = 10, maxSize = 1000, eps = 0) %>%
  lapply(arrange, pval)































