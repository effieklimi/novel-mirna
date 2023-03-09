library(tidyverse)
library(dplyr)
library("GSA")
library(fgsea)
library(plyr)
library(grid)

setwd("/Users/effieklimi/Documents/novel-mirna/")

########## pathways ##########
bp <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021"))
ke <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"))
re <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022"))
all_paths <- c(bp, ke, re)
names(all_paths) <- tolower(names(all_paths))
##############################

############## files ##############
# Open files containg gene and DE info per miRNA -> put them all inside a list:
deseqFiles <-
  list.files("/deseqResults/diffExprTables",
    pattern = "*.csv", full.names = TRUE) %>%
  lapply(read_csv) %>%
  lapply(distinct, name, .keep_all = TRUE)

# Name tibbles based on their filenames (miRNA name)
names(deseqFiles) <-
  list.files("/deseqResults/diffExprTables",
    pattern = "*.csv", full.names = FALSE) %>%
  lapply(gsub, pattern = "vsmiRCTRL_noLFCthreshold.csv", replacement = "")
# Save list of dfs containing all DE genes (without duplicates):
saveRDS(deseqFiles, file = "/deseqResults/deseqResultsTibble.rds")
####################################

# Create gene lists with LFC and gene symbols as names, ranked from high LFC to low
geneLists <-
  lapply(deseqFiles, "[", , c(8, 4)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)
# Save list of DE genes per miRNA with LFC (it's ranked):
saveRDS(geneLists, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/GSEA/NoLFCthresh/geneLists.rds")

# Perform GSEA
fgseaResults <-
  map(geneLists, fgsea, pathways = all_paths, minSize = 20, maxSize = 1000, eps = 0) %>%
  saveRDS(fgseaResults, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/GSEA/NoLFCthresh/fgseaResults.rds")

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
  `colnames<-`(c(names(deseqFiles), "Database"))

#### Save table of results for heatmap: ####
write.csv(fgseaResultsSig, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/GSEA/NoLFCthresh/fgseaResultsSigNES.csv")
############################################