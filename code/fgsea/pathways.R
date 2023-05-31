library("tidyverse")
library("fgsea")
library("GSA")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
# --------------------------------------------------

# Filter out lowly expressed genes:
expressed <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")
# functions:
read.geneset <- function(path_to_gset)  {
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]] != ""])
  names(out) = bp$geneset.names
  return(out)
}

get_gset_names_by_category = function(cat, gsets){
  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}
 
# Get gene sets:
bp <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021"))
saveRDS(bp, file = "results/rds/pathways/pathways-all-bioprocess.rds")

ke <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"))
saveRDS(ke, file = "results/rds/pathways/pathways-all-kegg.rds")

re <- read.geneset(url("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022"))
saveRDS(re, file = "results/rds/pathways/pathways-all-reactome.rds")

paths <- c(bp, ke, re)
names(paths) <- tolower(names(paths))
pathsExpressed <- lapply(paths, function(x) x[x %in% expressed$name])
saveRDS(pathsExpressed, file = "results/rds/pathways/pathways-all-expressed.rds")

# get cell cycle-associated genesets
ccKeys = c("cell cycle", "proliferation", "cell division")
ccPaths = get_gset_names_by_category(ccKeys, tolower(names(pathsExpressed)))
ccGs = lapply(ccPaths, function(x) pathsExpressed[[x]][pathsExpressed[[x]] != ''])
names(ccGs) = ccPaths
pathsExpressed[['pathways']][['all_cc_associated']] = ccPaths
saveRDS(ccGs, file = "results/rds/pathways/pathways-cellcycle.rds")

# get migration-associated genesets
miKeys = c('migration', 'motility', 'adhesion', 'locomotion', 'movement', 'actin', 'localization', "cytoskeleton", "cytoskeletal")
miPaths = get_gset_names_by_category(miKeys, tolower(names(pathsExpressed)))
miGs = lapply(miPaths, function(x) pathsExpressed[[x]][pathsExpressed[[x]]!=''])
names(miGs) = miPaths
pathsExpressed[['pathways']][['all_mi_associated']] = miPaths
saveRDS(miGs, file = "results/rds/pathways/pathways-motility.rds")

# get apoptosis-associated genesets
delKeys = c('apoptosis', 'apoptotic', 'senescence', 'senescent')
delPaths = get_gset_names_by_category(delKeys, tolower(names(pathsExpressed)))
delGs = lapply(delPaths, function(x) pathsExpressed[[x]][pathsExpressed[[x]]!=''])
names(delGs) = delPaths
pathsExpressed[['pathways']][['all_chol_associated']] = delPaths
saveRDS(delGs, file = "results/rds/pathways/pathway-deleterious.rds")
