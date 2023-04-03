library("tidyverse")
library("fgsea")
library("GSA")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
expressed <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")

ccGs <- readRDS("results/rds/pathways/pathways-cellcycle.rds")
miGs <- readRDS("results/rds/pathways/pathways-motility.rds")
delGs <- readRDS("results/rds/pathways/pathway-deleterious.rds")

shrinkResults <-
  readRDS("results/rds/vsmc-deseq2.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
  lapply(filter, name %in% expressed$name) # remove genes not expressed

geneLists <-
  lapply(shrinkResults, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)

ccFgsea <-
  map(geneLists, fgsea, pathways = ccGs, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults, file = "results/rds/fgsea/vsmc-fgsea-cellcycle.rds")

miFgsea <-
  map(geneLists, fgsea, pathways = miGs, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults, file = "results/rds/fgsea/vsmc-fgsea-motility.rds")

delFgsea <-
  map(geneLists, fgsea, pathways = delGs, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults, file = "results/rds/fgsea/vsmc-fgsea-deleterious.rds")


# Filter for p value, keep the Normalised Enrichment Scores and format table for the heatmap:
ccFgseaSig <-
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(ccFgsea, filter, padj < 0.05) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

miFgseaSig <-
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(miFgsea, filter, padj < 0.05) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

delFgseaSig <-
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(delFgsea, filter, padj < 0.05) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

#### Save table of results for heatmap: ####
write.csv(ccFgseaSig, file = "results/tables/fgsea-vsmc-cellcycle-sigNES.csv")
write.csv(miFgseaSig, file = "results/tables/fgsea-vsmc-motility-sigNES.csv")
write.csv(delFgseaSig, file = "results/tables/fgsea-vsmc-deleterious-sigNES.csv")
############################################