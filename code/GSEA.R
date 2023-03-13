library(tidyverse)
library(dplyr)
library("GSA")
library(fgsea)
library(plyr)
library(grid)

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

########## pathways ##########
read.geneset <- function(path_to_gset){
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]] != ''])
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
  readRDS("results/tables/vsmc-deseq2-results.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) # remove duplicates

geneLists <-
  lapply(shrinkResults, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)
# Save list of DE genes per miRNA with LFC (it's ranked):
saveRDS(geneLists, file = "results/tables/geneLists.rds")

fgseaResults <-
  map(geneLists, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0)
saveRDS(fgseaResults, file = "results/tables/vsmc-fgsea.rds")

# Filter for p value, keep the Normalised Enrichment Scores and format table for the heatmap:
fgseaResultsSig <-
  lapply(fgseaResults, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(names(geneLists), "Database"))

#### Save table of results for heatmap: ####
write.csv(fgseaResultsSig, file = "results/tables/fgsea-results-sigNES.csv")
############################################