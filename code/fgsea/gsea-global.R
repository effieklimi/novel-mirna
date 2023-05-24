library(tidyverse)
library(dplyr)
library("GSA")
library(fgsea)
library(plyr)
library(grid)

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
expressed <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")

########## pathways ##########
read.geneset <- function(path_to_gset)  {
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]] != ""])
  names(out) = bp$geneset.names
  return(out)
}

bp <- readRDS("results/rds/pathways/pathways-all-bioprocess.rds") 
ke <- readRDS("results/rds/pathways/pathways-all-kegg.rds")
re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")
names(bp) %>% lapply(str_to_sentence) %>% lapply(gsub, pattern = "Dna", replacement = "DNA", fixed = TRUE) -> names(bp)

pathsExpressed <- c(bp, ke, re)
names(pathsExpressed) <- tolower(names(pathsExpressed))

goReactome <- c(bp, re)
##############################

shrinkResults <-
  readRDS("results/rds/vsmc-deseq2.rds") %>%
  lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
  lapply(filter, name %in% expressed$name) # remove genes not expressed

geneLists <-
  lapply(shrinkResults, "[", , c(7, 3)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = TRUE)
# Save list of DE genes per miRNA with LFC (it's ranked):
saveRDS(geneLists, file = "results/rds/genelists-vsmc.rds")

fgseaResults <-
  map(geneLists, fgsea, pathways = goReactome, minSize = 20, maxSize = 1000, eps = 0, nPermSimple = 100000)
# Discussion on nPerms: https://www.biostars.org/p/387492/

saveRDS(fgseaResults, file = "results/rds/vsmc-fgsea.rds")

# Filter for p value, keep the Normalised Enrichment Scores and format table for the heatmap:
fgseaResultsSig <-
  #lapply(fgseaResults, arrange, -pval) %>%
  lapply(fgseaResults, filter, padj < 0.01) %>%
  lapply("[", , c("pathway", "NES")) %>%
  join_all(by = "pathway", type = "left") %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c(miRNAnames))

#### Save table of results for heatmap: ####
write.csv(fgseaResultsSig, file = "results/tables/fgsea-vsmc-sigNES.csv")
############################################