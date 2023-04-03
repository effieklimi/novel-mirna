library(tidyverse)
library(dplyr)
library(GSA)
library(fgsea)
library(plyr)
library(grid)


options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

########## pathways ##########
read.geneset <- function(path_to_gset)  {
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]] != ""])
  names(out) = bp$geneset.names
  return(out)
}

bp <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021"))
ke <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"))
re <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022"))
all_paths <- c(bp, ke, re)
names(all_paths) <- tolower(names(all_paths))
##############################

vsmc50Top2 <-
  readRDS("results/rds/p01/targets50top2-vsmc-p01.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)

geneLists50Top2 <-
vsmc50Top2 %>%
lapply(tibble::deframe) %>%
lapply(sort, decreasing = TRUE)

fgseaResults50Top2 <-
  map(geneLists50Top2, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)

# Filter for p value, keep the Normalised Enrichment Scores and format table for the heatmap:
fgseaResultsSig <- 
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(fgseaResults, filter, padj < 0.001) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

#### Save table of results for heatmap: ####
write.csv(fgseaResultsSig, file = "results/tables/fgseaResultsSig50Top2NES-multimiR-p01-vsmc.csv")
############################################











































options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

read.geneset <- function(path_to_gset) {
  bp <- GSA.read.gmt(path_to_gset)
  out <- bp$genesets
  out <- lapply(1:length(out), function(x) out[[x]][out[[x]] != ""])
  names(out) <- bp$geneset.names
  return(out)
}


bp <- read.geneset("databases/GO_Biological_Process_2021.txt")
ke <- read.geneset("databases/KEGG_2021_Human.txt")
re <- read.geneset("databases/Reactome_2022.txt")
all_paths <- c(bp, ke, re)
names(all_paths) <- tolower(names(all_paths))


vsmc50Top2 <-
  readRDS("results/rds/p01/targets50top2-vsmc-p01.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)

geneLists50Top2 <-
  vsmc50Top2 %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)

fgseaResults50Top2 <-
  map(geneLists50Top2, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)

fgseaResults50Top2 <-
  lapply(fgseaResults50Top2, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  `colnames<-`(c("pathway", miRNAnames)) %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway")


fgseaResults50Top2 <-
  lapply(arrange, -pval)
  lapply(fgseaResults50Top2, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>% 
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(names(geneLists)))
  

write.csv(fgseaResultsSig50Top2, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/fgseaResultsSig50Top2NES-multimiR-p01-vsmc.csv")












targets50Top2 <-
readRDS("/Users/effieklimi/Documents/novel-mirna/results/tables/targets50Top2.rds") %>%
  lapply(distinct, name, .keep_all = TRUE)
targets50 <-
readRDS("/Users/effieklimi/Documents/novel-mirna/results/tables/targets50.rds")  %>%
  lapply(distinct, name, .keep_all = TRUE)
targets100 <-
readRDS("/Users/effieklimi/Documents/novel-mirna/results/tables/targets100.rds")  %>%
  lapply(distinct, name, .keep_all = TRUE)




geneLists50Top2 <-
  lapply(targets50Top2, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)

geneLists50 <-
  lapply(targets50, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)

geneLists100 <-
  lapply(targets100, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)




fgseaResults50Top2 <-
  map(geneLists50Top2, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults50Top2, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResults50Top2.rds")

fgseaResults50 <-
  map(geneLists50, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults50, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResults50.rds")

fgseaResults100 <-
  map(geneLists100, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults100, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResults100.rds")


fgseaResultsSig50Top2 <-
  lapply(fgseaResults50Top2, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  `colnames<-`(c("pathway", names(deseqFiles))) %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway")

fgseaResultsSig50 <-
  lapply(fgseaResults50, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  `colnames<-`(c("pathway", names(deseqFiles))) %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway")

fgseaResultsSig100 <-
  lapply(fgseaResults100, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  `colnames<-`(c("pathway", names(deseqFiles))) %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway")


write.csv(fgseaResultsSig50Top2, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/fgseaResultsSig50Top2NES.csv")
write.csv(fgseaResultsSig50, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/ffgseaResultsSig50NES.csv")
write.csv(fgseaResultsSig100, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/ffgseaResultsSig100NES.csv")