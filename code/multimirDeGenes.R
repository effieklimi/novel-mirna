library(tidyverse)
library(dplyr)
library(multiMiR)

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

ensIdSplit <- function(x) { strsplit(x, split = ".", fixed = TRUE)[[1]][1] }
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
deseqFiles <- readRDS("results/tables/vsmc-deseq2-results.rds")
deGenesList <- 
  lapply(deseqFiles, "[", , c(1, 3, 7)) %>%
  map(~ mutate(.x, EnsID = ensIdSplit(EnsID)))



# MultimiR for all 7 miRNAs, top 50% of at least 2 tools
multimirPredicted50Freq2 <- map(
  miRNAnames, 
  ~ get_multimir(
    org = 'hsa', 
    table = 'predicted', 
    mirna = ., 
    predicted.cutoff.type = "p", 
    predicted.cutoff = 50, 
    predicted.site = "all"
  )
) %>% 
  lapply(function(x) x@data) %>%
  lapply("[", , 4) %>%
  lapply(unlist) %>%
  lapply(table) %>%
  lapply(as.data.frame) %>%
  lapply(filter, Freq >= 2)


# MultimiR for all 7 miRNAs, top 50% of any
multimirPredicted50 <- map(
  miRNAnames, 
  ~ get_multimir(
    org = 'hsa', 
    table = 'predicted', 
    mirna = ., 
    predicted.cutoff.type = "p", 
    predicted.cutoff = 50, 
    predicted.site = "all"
    )
  ) %>% 
  lapply(function(x) x@data) %>%
  lapply("[", , c(6, 4)) 




# MultimiR for all 7 miRNAs without using prediction score cutoffs
multimirPredicted100 <- map(
  miRNAnames, 
  ~ get_multimir(
    org = 'hsa', 
    table = 'predicted', 
    mirna = .,    
    predicted.cutoff.type = "n", 
    predicted.cutoff = 100000000000, 
    predicted.site = "all"
    )
  ) %>% 
  lapply(function(x) x@data) %>%
  lapply("[", , c(6, 4))
  

targets50Top2 <- 
  mapply(function(x,y) x[which(as.list(x$name) %in% y$Var1),], deGenesList, multimirPredicted50Freq2, SIMPLIFY = F) %>%
  lapply(filter, log2FoldChange < 0)

targets50 <- 
  mapply(function(x,y) x[which(as.list(x$name) %in% as.list(y$target_symbol)),], deGenesList, multimirPredicted50, SIMPLIFY = F) %>%
  lapply(filter, log2FoldChange < 0)

targets100 <- 
  mapply(function(x,y) x[which(as.list(x$name) %in% as.list(y$target_symbol)),], deGenesList, multimirPredicted100, SIMPLIFY = F) %>%
  lapply(filter, log2FoldChange < 0)


names(targets50Top2) <- names(deseqFiles)
names(targets50) <- names(deseqFiles)
names(targets100) <- names(deseqFiles)

saveRDS(targets50Top2, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/targets50Top2.rds")
saveRDS(targets50, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/targets50.rds")
saveRDS(targets100, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/targets100.rds")



