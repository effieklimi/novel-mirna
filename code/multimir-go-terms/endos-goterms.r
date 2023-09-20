library("clusterProfiler")
library("tidyverse")
library("org.Hs.eg.db")


endosDeseq <-
  readRDS("results/rds/p01/endos-deseq2-p01.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) # remove duplicates

endosGenes <-
  lapply(endosDeseq, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE) %>%
  lapply(tibble::enframe)

  
targets50top2 <- readRDS("/Users/effieklimi/Documents/novel-mirna/results/rds/p01/targets50Top2-endos-p01.rds")
FPKM <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/endosFpkm.csv", header = TRUE)

FPKM <- distinct(FPKM, name, .keep_all = TRUE)
FPKMrownmeans <- rowMeans(FPKM[,c(8:40)], na.rm = TRUE)
FPKM$rowMeans <- FPKMrownmeans
FPKM <- filter(FPKM, rowMeans > 2)
# For the pathway analysis we took the uppermost quartile wrt logfoldchange
geneLists50Top2BP <-
  lapply(targets50Top2, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  #lapply(tibble::enframe) %>%
  #lapply(filter, value < quantile(value, .25)) %>%
  #lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO,
      universe      = FPKM$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1) %>%
  lapply(function(x) x@result)

pdf(file = "results/figures/endos-mir323-50Top2BP-p01.pdf", width = 8, height = 4.5)
geneLists50Top2BP$condition_mir323_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) +
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2, 4, 6, 8, 10)) +
  scale_x_discrete(labels = label_wrap(60)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-323-3p: GO Biological Process")
dev.off()

pdf(file = "results/figures/endos-mir892-50Top2BP-p01.pdf", width = 8, height = 4.5)
geneLists50Top2BP$condition_mir892b_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) +
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2, 4, 6, 8, 10)) +
  scale_x_discrete(labels = label_wrap(65)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-892b: GO Biological Process")
dev.off()