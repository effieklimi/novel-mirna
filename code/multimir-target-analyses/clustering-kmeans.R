library("ComplexHeatmap")
library("tidyverse")
library("org.Hs.eg.db")
library("enrichR")
library("clusterProfiler")
library(dendextend)
library(plyr)


options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
ensIdSplit <- function(x) { strsplit(x, split = ".", fixed = TRUE)[[1]][1] }

# naming lists
miRnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
colNames <- c("miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")
#



# 3. hclust -> lfc filt post hoc -> go enrichment


countFpkm <- read.csv(
  "results/tables/vsmcFpkm.csv",
  header = TRUE
) %>% as_tibble()

deseq <- # "deseq" to be used for clustering using DE only (no multimir)
  readRDS("results/rds/vsmc-deseq2.rds") %>%
  lapply("[", , c(1, 3, 7)) %>%
  map(~ mutate(.x, EnsID = ensIdSplit(EnsID))) %>%
  lapply(as_tibble) %>%
  lapply(filter, name %in% countFpkm$name) %>%
  lapply(filter, log2FoldChange < -log2(1.5))
names(deseq) <- miRnames
deseq <- deseq[c(2, 3, 5, 7, 1, 4, 6)]
targets <- do.call(rbind, deseq)

# loading expression data and multimir data
multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(filter, log2FoldChange < -log2(1.5))
names(multimir) <- miRnames
multimir <- multimir[c(2, 3, 5, 7, 1, 4, 6)]
targets <- do.call(rbind, multimir)

genes <- join(targets[, 3], countFpkm[, c(2, 17:40)], by = "name") %>% na.omit %>% distinct(name, .keep_all = TRUE)
matrix <- genes[, 2:25]
rownames(matrix) <- genes$name
matrix <- as.matrix(matrix) %>% +1 %>% log2()
mean <- rowMeans(matrix)
sd <- apply(matrix, 1, sd)
matrix <- (matrix - mean) / sd 
colnames(matrix) <- colNames
#

# clustering kmeans
matrix <- matrix %>% na.omit()
clusterNo <- 6
kmRes <- kmeans(matrix, clusterNo, nstart = 25, iter.max = 30)

first <-  names(kmRes$cluster[kmRes$cluster == 1])
second <- names(kmRes$cluster[kmRes$cluster == 2])
third <- names(kmRes$cluster[kmRes$cluster == 3])
fourth <- names(kmRes$cluster[kmRes$cluster == 4])
fifth <- names(kmRes$cluster[kmRes$cluster == 5])
sixth <- names(kmRes$cluster[kmRes$cluster == 6])

#clusters <- list(first, second, third, fourth)
clusters <- list(first, second, third, fourth, fifth, sixth)
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
pdf("results/figures/deseq-heatmap-clusters-kmeans-6.pdf", width = 7, height = 10)
ComplexHeatmap::pheatmap(as.matrix(clusterMatrix),
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 8,
    show_rownames = FALSE,
    show_colnames = TRUE,
    row_split = rep(clustersNo$name, clustersNo$geneNo),
    legend = TRUE,
    annotation_legend = TRUE,
    colorRampPalette(c("#003ba1", "#ffffff", "#950046"))(500)
)
dev.off()
#

# no padj filtering during the analysis, so I can explore based on ranking only for starters
# If needed, I can filter based on padj after
geneOntology <-
  clusters %>%
  map(enrichGO,
      universe      = countFpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.1,
      qvalueCutoff  = 1,
      minGSSize     = 20)

geneOntologyTable <- lapply(geneOntology, function(x) x@result)
saveRDS(geneOntologyTable, file = "results/rds/vsmc-deseq-gobp-kmeans-6-lfc-filtered.rds")

bpSimplify <- lapply(geneOntology, simplify, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
bpSimpTable <- lapply(bpSimplify, function(x) x@result)
saveRDS(bpSimplify, file = "results/rds/vsmc-deseq-gobp-kmeans-6-lfc-filtered-simplified.rds")
#



# GO term dot plots
pdf(file = "results/figures/vsmc-deseq-gobp-cluster1-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[1]], showCategory = 20) +
ggplot2::ggtitle("Cluster 1") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()


pdf(file = "results/figures/vsmc-deseq-gobp-cluster2-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[2]], showCategory = 20) +
ggplot2::ggtitle("Cluster 2") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    legend.key.size = unit(1, 'cm'),
    title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()


pdf(file = "results/figures/vsmc-deseq-gobp-cluster3-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[3]], showCategory = 20) +
ggplot2::ggtitle("Cluster 3") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    legend.key.size = unit(1, 'cm'),
    title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()


pdf(file = "results/figures/vsmc-deseq-gobp-cluster4-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[4]], showCategory = 20) +
ggplot2::ggtitle("Cluster 4") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    legend.key.size = unit(1, 'cm'),
    title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()

pdf(file = "results/figures/vsmc-deseq-gobp-cluster5-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[5]], showCategory = 20) +
ggplot2::ggtitle("Cluster 5") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    legend.key.size = unit(1, 'cm'),
    title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()

pdf(file = "results/figures/vsmc-deseq-gobp-cluster6-kmeans6.pdf", width = 9, height = 12)
dotplot(bpSimplify[[6]], showCategory = 20) +
ggplot2::ggtitle("Cluster 6") +
  ggplot2::xlab("Enrichment") +
  ggplot2::scale_color_gradient(low = "#ffc2df", high = "#7c0046") +
  theme(
    legend.key.size = unit(1, 'cm'),
    title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 10, face = "bold"))
dev.off()

# Cluster target percentage

## ~~~ miR-323a-3p
mir323Perc <- data.frame(
  "miRNA" = c("miR-323a-3p", "miR-323a-3p", "miR-323a-3p", "miR-323a-3p"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[1]], name %in% clusters[[1]])) / nrow(multimir[[1]])) * 100,
      (nrow(filter(multimir[[1]], name %in% clusters[[2]])) / nrow(multimir[[1]])) * 100,
      (nrow(filter(multimir[[1]], name %in% clusters[[3]])) / nrow(multimir[[1]])) * 100,
      (nrow(filter(multimir[[1]], name %in% clusters[[4]])) / nrow(multimir[[1]])) * 100)
  )

## ~~~ miR-449b-5p
mir449Perc <- data.frame(
  "miRNA" = c("miR-449b-5p", "miR-449b-5p", "miR-449b-5p", "miR-449b-5p"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[2]], name %in% clusters[[1]])) / nrow(multimir[[2]])) * 100,
      (nrow(filter(multimir[[2]], name %in% clusters[[2]])) / nrow(multimir[[2]])) * 100,
      (nrow(filter(multimir[[2]], name %in% clusters[[3]])) / nrow(multimir[[2]])) * 100,
      (nrow(filter(multimir[[2]], name %in% clusters[[4]])) / nrow(multimir[[2]])) * 100)
  )

## ~~~ miR-491-3p
mir491Perc <- data.frame(
  "miRNA" = c("miR-491-3p", "miR-491-3p", "miR-491-3p", "miR-491-3p"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[3]], name %in% clusters[[1]])) / nrow(multimir[[3]])) * 100,
      (nrow(filter(multimir[[3]], name %in% clusters[[2]])) / nrow(multimir[[3]])) * 100,
      (nrow(filter(multimir[[3]], name %in% clusters[[3]])) / nrow(multimir[[3]])) * 100,
      (nrow(filter(multimir[[3]], name %in% clusters[[4]])) / nrow(multimir[[3]])) * 100)
  )

## ~~~ miR-892b
mir892bPerc <- data.frame(
  "miRNA" = c("miR-892b", "miR-892b", "miR-892b", "miR-892b"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[4]], name %in% clusters[[1]])) / nrow(multimir[[4]])) * 100,
      (nrow(filter(multimir[[4]], name %in% clusters[[2]])) / nrow(multimir[[4]])) * 100,
      (nrow(filter(multimir[[4]], name %in% clusters[[3]])) / nrow(multimir[[4]])) * 100,
      (nrow(filter(multimir[[4]], name %in% clusters[[4]])) / nrow(multimir[[4]])) * 100)
  )

## ~~~ miR-1827
mir1827Perc <- data.frame(
  "miRNA" = c("miR-1827", "miR-1827", "miR-1827", "miR-1827"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[5]], name %in% clusters[[1]])) / nrow(multimir[[5]])) * 100,
      (nrow(filter(multimir[[5]], name %in% clusters[[2]])) / nrow(multimir[[5]])) * 100,
      (nrow(filter(multimir[[5]], name %in% clusters[[3]])) / nrow(multimir[[5]])) * 100,
      (nrow(filter(multimir[[5]], name %in% clusters[[4]])) / nrow(multimir[[5]])) * 100)
  )

## ~~~ miR-4774-3p
mir4774Perc <- data.frame(
  "miRNA" = c("miR-4774-3p", "miR-4774-3p", "miR-4774-3p", "miR-4774-3p"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[6]], name %in% clusters[[1]])) / nrow(multimir[[6]])) * 100,
      (nrow(filter(multimir[[6]], name %in% clusters[[2]])) / nrow(multimir[[6]])) * 100,
      (nrow(filter(multimir[[6]], name %in% clusters[[3]])) / nrow(multimir[[6]])) * 100,
      (nrow(filter(multimir[[6]], name %in% clusters[[4]])) / nrow(multimir[[6]])) * 100)
  )

## ~~~ miR-5681b
mir5681bPerc <- data.frame(
  "miRNA" = c("miR-5681b", "miR-5681b", "miR-5681b", "miR-5681b"),
  "Clusters" = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  "Percentage" = c(
      (nrow(filter(multimir[[7]], name %in% clusters[[1]])) / nrow(multimir[[7]])) * 100,
      (nrow(filter(multimir[[7]], name %in% clusters[[2]])) / nrow(multimir[[7]])) * 100,
      (nrow(filter(multimir[[7]], name %in% clusters[[3]])) / nrow(multimir[[7]])) * 100,
      (nrow(filter(multimir[[7]], name %in% clusters[[4]])) / nrow(multimir[[7]])) * 100)
  )

percentage <- rbind(mir323Perc, mir449Perc, mir491Perc, mir892bPerc, mir1827Perc, mir4774Perc, mir5681bPerc)
pdf("results/figures/multimir-cluster-percetages.pdf", width = 8, height = 4)
ggplot(percentage, aes(fill=Clusters, y=Percentage, x = factor(miRNA, level = c('miR-323a-3p', 'miR-449b-5p', 'miR-491-3p', 'miR-892b', 'miR-1827', 'miR-4774-3p', 'miR-5681b')))) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = c("#6288c5", "#adbbe4", "#e989c0", "#e05694")) +
    theme_bw() +
    xlab(' ') +
    ylab("Percentage of targets in each cluster") +
    scale_y_continuous(labels = scales::percent)
dev.off()