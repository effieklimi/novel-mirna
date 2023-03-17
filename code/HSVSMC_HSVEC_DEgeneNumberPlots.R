library(reshape2)
library(viridis)
library(ReactomePA)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)




######### FPKM files #########
fpkmEC <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/endosFpkm.csv", header = TRUE)
fpkmSMC <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcFpkm.csv", header = TRUE)

fpkmEC <- distinct(fpkmEC, name, .keep_all = TRUE)
fpkmECrownmeans <- rowMeans(fpkmEC[,c(8:37)], na.rm = TRUE)
fpkmEC$rowMeans <- fpkmECrownmeans
fpkmEC <- filter(fpkmEC, rowMeans > 2)

fpkmSMC <- distinct(fpkmSMC, name, .keep_all = TRUE)
fpkmSMCrownmeans <- rowMeans(fpkmSMC[,c(8:40)], na.rm = TRUE)
fpkmSMC$rowMeans <- fpkmSMCrownmeans
fpkmSMC <- filter(fpkmSMC, rowMeans > 2)

# Overlap of the total transcriptomes:
cols <- brewer.pal(2, "Pastel2")
venn.diagram(
  x = list(fpkmSMC$ENSEMBL, fpkmEC$ENSEMBL),
  category.names = c("SMCs" , "ECs"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#CBD5E8", "#FDCDAC"),
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans")

############# DE files #############
# Open files containg gene and DE info per miRNA -> put them all inside a list:

endosDeseq <-
  readRDS("results/rds/p01/endos-deseq2-p01.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) # remove duplicates

vsmcDeseq <-
  readRDS("results/rds/p01/vsmc-deseq2-p01.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) # remove duplicates

endosGenes <-
  lapply(endosDeseq, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE) %>%
  lapply(tibble::enframe)

vsmcGenes <-
  lapply(vsmcDeseq, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE) %>%
  lapply(tibble::enframe)


#### Data frames for stakced bar plots of DE genes (up and down)
###### 
deGenesSMC <-
  data.frame(
    upreg = unlist(lapply(vsmcGenes, subset, value > 0) %>% lapply(nrow)),
    downreg = -unlist(lapply(vsmcGenes, subset, value < 0) %>% lapply(nrow))
  )


deGenesEC <-
  data.frame(
    upreg = unlist(lapply(endosGenes, subset, value > 0) %>% lapply(nrow)),
    downreg = -unlist(lapply(endosGenes, subset, value < 0) %>% lapply(nrow))
  )


rownames(deGenesEC) <- c(
  "miR-1827",
  "miR-323a-3p",
  "miR-449b-5p",
  "miR-4774-3p",
  "miR-491-3p",
  "miR-5681b",
  "miR-892b"
)

rownames(deGenesSMC) <- c(
  "miR-1827",
  "miR-323a-3p",
  "miR-449b-5p",
  "miR-4774-3p",
  "miR-491-3p",
  "miR-5681b",
  "miR-892b"
)
###############################


EC <- t(deGenesEC)[, c(
  "miR-5681b",
  "miR-4774-3p",
  "miR-1827",
  "miR-892b",
  "miR-491-3p",
  "miR-449b-5p",
  "miR-323a-3p"
)] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
  melt()
# with control
ggplot(EC, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#eaa072", "#72bcea")) +
  theme_minimal() +
  xlab("Conditions") +
  ylab("Number of Genes") +
  coord_flip()
# without control
ggplot(EC, aes(y = value, x = sample, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#eaa072", "#72bcea")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of DE Genes") +
  coord_flip()


# xxxxxxxxxxxxxxxxxxxx #

SMC <-
  t(deGenesSMC)[, c(
    "miR-5681b",
    "miR-4774-3p",
    "miR-1827",
    "miR-892b",
    "miR-491-3p",
    "miR-449b-5p",
    "miR-323a-3p"
)] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  melt()
# stacked barplot
ggplot(SMC, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#eaa072", "#72bcea")) +
  theme_minimal() +
  xlab("Conditions") +
  ylab("Number of Genes") +
  coord_flip()
# perc barplot
ggplot(SMC, aes(y = value, x = sample, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#eaa072", "#72bcea")) +
  theme_minimal() +
  xlab("Conditions") +
  ylab("Number of Genes") +
  coord_flip()




#### Controls: EC
ECquiescGoTerms <- 
  enrichGO(
    gene = deseqFilesEC$FBS02vsFBS10$name,
    universe      = fpkmEC$name,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",    
    pAdjustMethod = "BH",
    pvalueCutoff  = 1)

ECquiescGoTerms@result[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=1, width=.6) + 
  coord_flip() +
  theme_minimal() +
  xlab(" ") +
  labs(title = "GO TERMS (BP)",
       subtitle = "DE genes in response to 10%FBS",
       caption = "Top 20 most significant GO terms") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0)
  )

fgseaResultsEC <- fgsea(
  geneListsEC$FBS02vsFBS10, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0
  ) # Nothing sig

##############################################################
##############################################################
##############################################################


#### Controls: EC
##############################################################
##############################################################
SMCquiescGoTerms <- 
  enrichPathway(
    gene          = deseqFilesSMC$IPvsFBS02$name,
    universe      = fpkmSMC$name,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",    
    pAdjustMethod = "BH",
    pvalueCutoff  = 1)

SMCquiescGoTerms[c(1:20),] %>%
  arrange(-pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Description=factor(Description, levels=Description)) %>%   # This trick update the factor levels
  ggplot(aes(x=Description, y=-log(pvalue))) +
  geom_bar(stat="identity", fill="#f68060", alpha=1, width=.6) + 
  coord_flip() +
  theme_minimal() +
  labs(title = "GO TERMS (BP)",
       subtitle = "DE genes in response to IL-1A/PDGF-BB",
       caption = "Top 20 most significant GO terms") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = -1)
  )

fgseaResultsSMC <- fgsea(
  geneListsSMC$IPvsFBS02, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0
)
