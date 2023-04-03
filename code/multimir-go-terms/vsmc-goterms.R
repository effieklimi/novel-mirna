library("clusterProfiler")
library("tidyverse")
library("purrr")
library("org.Hs.eg.db")
library("scales")
library("reshape2")
library("ggbreak")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
targets <- readRDS("results/rds/vsmc-multimir.rds")

bpOrder <- c(targets[2], targets[3], targets[5], targets[7], targets[1], targets[4], targets[6])


fpkm <- read.csv(
  "results/tables/vsmcExpressed.csv",
  header = TRUE
) %>% as_tibble()

# For the pathway analysis we took the uppermost quartile wrt logfoldchange
bpEnrichment <-
  lapply(bpOrder, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .25)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO,
      universe      = fpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BY",
      pvalueCutoff  = 1) %>%
  lapply(function(x) x@result)

saveRDS(bpEnrichment, file = "results/rds/vsmc-multimir-gobp.rds")

bpOrder <- lapply(bpEnrichment, function(x) rbind(head(x, 10)))
bpOrder <- 
  lapply(bpOrder, "[", c(2, 5)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) 

results <- data.frame(
  mir323 = c(bpOrder[[1]]$value, rep(NA, 60)),
  mir449b = c(rep(NA, 10), bpOrder[[2]]$value, rep(NA, 50)),
  mir491 = c(rep(NA, 20), bpOrder[[3]]$value, rep(NA, 40)),
  mir892b = c(rep(NA, 30), bpOrder[[4]]$value, rep(NA, 30)),
  mir1827 = c(rep(NA, 40), bpOrder[[5]]$value, rep(NA, 20)),
  mir4774 = c(rep(NA, 50), bpOrder[[6]]$value, rep(NA, 10)),
  mir5681b = c(rep(NA, 60), bpOrder[[7]]$value)
)
results <- results[order(nrow(results):1),]

names = c(rev(bpOrder[[7]]$name), 
          rev(bpOrder[[6]]$name), 
          rev(bpOrder[[5]]$name), 
          rev(bpOrder[[4]]$name), 
          rev(bpOrder[[3]]$name), 
          rev(bpOrder[[2]]$name), 
          rev(bpOrder[[1]]$name))
names <- make.names(names, unique = TRUE)
names <- gsub(".", " ", names, fixed=TRUE)

results$names <- names
colnames(results) <- c("miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b", "name")






write.csv(results, file = "results/tables/vsmc-multimir-gobp.csv")

resultsMelt <- melt(results, id.vars = "name")
resultsMelt$name <- factor(names, levels = names)

pdf(file = "results/figures/vsmc-multimir-gobp.pdf", width = 12, height = 10)
ggplot(data=resultsMelt) + 
  aes(x=variable, y=name, size=-log(value)) + 
  geom_point(alpha=0.5, col='#b7006e') +
  theme_light() 
dev.off()








































geneLists50Top2MF <-
  lapply(targets, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .25)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO, 
      universe      = fpkm$name,
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