library("clusterProfiler")
library("tidyverse")
library("purrr")
library("org.Hs.eg.db")
library("scales")

targets50Top2 <- readRDS("/Users/effieklimi/Documents/novel-mirna/results/rds/p01/targets50Top2-vsmc-p01.rds")

FPKM <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcFpkm.csv", header = TRUE)
FPKM <- distinct(FPKM, name, .keep_all = TRUE)
FPKMrownmeans <- rowMeans(FPKM[,c(8:40)], na.rm = TRUE)
FPKM$rowMeans <- FPKMrownmeans
FPKM <- filter(FPKM, rowMeans > 2)


nrowfun   <- function(vector, threshold){
  num_to_thres <- floor(threshold*0.01*length(vector))
  l = length (vector)
  score = c(rep("H",num_to_thres),rep("L",l-num_to_thres))
  return(score)
}



# For the pathway analysis we took the uppermost quartile wrt logfoldchange
geneLists50Top2BP <-
  lapply(targets50Top2, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .25)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO,
      universe      = FPKM$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1) %>%
  lapply(function(x) x@result)

geneLists50Top2MF <-
  lapply(targets50Top2, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .25)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO, 
      universe      = FPKM$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1) %>%
  lapply(function(x) x@result)






#BP
geneLists50Top2BP$condition_mir1827_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(65)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-1827: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir1827-50Top2BP-p01.pdf", width = 8, height = 5.5)


geneLists50Top2BP$condition_mir323_vs_mirctrl[c(1:10), ] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description = factor(Description, levels = Description)) %>%   # This trick update the factor levels
  ggplot(aes(x = Description, y = -log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(65)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-323-3p: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir323-50Top2BP-p01.pdf", width = 8, height = 5.5)


geneLists50Top2BP$condition_mir449b_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0,4,8,12, 16)) +
  scale_x_discrete(labels = label_wrap(50)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-449b-5p: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir449b-50Top2BP-p01.pdf", width = 7, height = 5.5)


geneLists50Top2BP$condition_mir4774_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(70)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-4774-3p: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir4774-50Top2BP-p01.pdf", width = 8, height = 5.5)


geneLists50Top2BP$condition_mir491_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(40)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-491-3p: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir491-50Top2BP-p01.pdf", width = 6, height = 5.5)


geneLists50Top2BP$condition_mir5681b_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by pvalue This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
 geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(50)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-5681b: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir5681b-50Top2BP-p01.pdf", width = 7, height = 5.5)


geneLists50Top2BP$condition_mir892b_vs_mirctrl[c(1:10),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", color = "black", fill="#5c96cc", alpha = 1, width=.7) + 
  coord_flip() +
  scale_y_discrete(expand = expansion(mult = c(0, .1)), limits = c(0, 2,4,6,8, 10)) +
  scale_x_discrete(labels = label_wrap(55)) +
  theme_minimal() +
  theme(axis.text = element_text(face="bold")) +
  theme(axis.line.x.bottom = element_line(size = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("miR-892b: GO Biological Process")
ggsave("/Users/effieklimi/Documents/novel-mirna/results/figures/mir892b-50Top2BP-p01.pdf", width = 7.5, height = 5.5)







#MP
geneLists50Top2MF$condition_mir1827_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-1827 Molecular Function")

geneLists50Top2MF$condition_mir323_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-323-3p Molecular Function")

geneLists50Top2MF$condition_mir449b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-449b-5p Molecular Function")

geneLists50Top2MF$condition_mir4774_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-4774-3p Molecular Function")

geneLists50Top2MF$condition_mir491_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-R491-3p Molecular Function")

geneLists50Top2MF$condition_mir5681b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by pvalue This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-5681b Molecular Function")

geneLists50Top2MF$condition_mir892b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-892b Molecular Function")