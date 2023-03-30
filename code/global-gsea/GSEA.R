library(tidyverse)
library(dplyr)
library("GSA")
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

shrinkResults <- 
  readRDS("results/rds/p01/vsmc-deseq2-p01.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) # remove duplicates

geneLists <-
  lapply(shrinkResults, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)
# Save list of DE genes per miRNA with LFC (it's ranked):
saveRDS(geneLists, file = "results/rds/p01/genelists-vsmc-p01.rds")

fgseaResults <-
  map(geneLists, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults, file = "results/rds/p01/vsmc-fgsea-p01.rds")

# Filter for p value, keep the Normalised Enrichment Scores and format table for the heatmap:
fgseaResultsSig <- 
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(fgseaResults, filter, padj < 0.001) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>% 
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

#### Save table of results for heatmap: ####
write.csv(fgseaResultsSig, file = "results/tables/fgsea-vsmc-sigNES-p01.csv")
############################################