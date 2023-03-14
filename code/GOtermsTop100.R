targets100 <- readRDS("/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/targets50Top2.rds")

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
geneLists100BP <-
  lapply(targets100, "[", , c(3, 2)) %>%
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

geneLists100MF <-
  lapply(targets100, "[", , c(3, 2)) %>%
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
geneLists100BP$condition_mir1827_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-1827 Biological Process")

geneLists100BP$condition_mir323_vs_mirctrl[c(1:20), ] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description = factor(Description, levels = Description)) %>%   # This trick update the factor levels
  ggplot(aes(x = Description, y = -log(pvalue))) +
  geom_bar(stat = "identity", fill = "#f68060", alpha = .6, width = .7) +
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-323-3p Biological Process")

geneLists100BP$condition_mir449b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-449b-5p Biological Process")

geneLists100BP$condition_mir4774_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-4774-3p Biological Process")

geneLists100BP$condition_mir491_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-R491-3p Biological Process")

geneLists100BP$condition_mir5681b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by pvalue This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-5681b Biological Process")

geneLists100BP$condition_mir892b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-892b Biological Process")






#MP
geneLists100MF$condition_mir1827_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-1827 Molecular Function")

geneLists100MF$condition_mir323_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-323-3p Molecular Function")

geneLists100MF$condition_mir449b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-449b-5p Molecular Function")

geneLists100MF$condition_mir4774_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal() +
  ggtitle("miR-4774-3p Molecular Function")

geneLists100MF$condition_mir491_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-R491-3p Molecular Function")

geneLists100MF$condition_mir5681b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by pvalue This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-5681b Molecular Function")

geneLists100MF$condition_mir892b_vs_mirctrl[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.7) + 
  coord_flip() +
  theme_minimal()+
  ggtitle("miR-892b Molecular Function")