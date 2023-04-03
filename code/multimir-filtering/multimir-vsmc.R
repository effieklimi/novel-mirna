library(tidyverse)
library(dplyr)
library(multiMiR)

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

ensIdSplit <- function(x) { strsplit(x, split = ".", fixed = TRUE)[[1]][1] }
mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

fpkm <- read.csv(
  "results/tables/vsmcFpkm.csv",
  header = TRUE
) %>% as_tibble()

fpkm <- distinct(fpkm, name, .keep_all = TRUE)
fpkm <- filter(
  fpkm, quiesmean > 2 |
  ipmean > 2 |
  mockmean > 2 |
  mirctrlmean > 2 |
  mir323amean > 2 |
  mir449bmean > 2 |
  mir491mean > 2 |
  mir892mean > 2 |
  mir1827mean > 2 |
  mir4774mean > 2 |
  mir5681bmean > 2
)

write.csv(fpkm, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")

deseq <- readRDS("results/rds/vsmc-deseq2.rds")
genes <- 
  lapply(deseq, "[", , c(1, 3, 7)) %>%
  map(~ mutate(.x, EnsID = ensIdSplit(EnsID))) %>%
  lapply(as_tibble) %>%
  lapply(filter, name %in% fpkm$name)
names(genes) <- mirNames

# MultimiR for all 7 miRNAs, top 50% of at least 2 tools
multimirPredicted50Freq2 <- map(
  mirNames, 
  ~ get_multimir(
    org = "hsa",
    table = "predicted",
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

targets50Top2 <- 
  mapply(function(x, y) x[which(as.list(x$name) %in% y$Var1),], genes, multimirPredicted50Freq2, SIMPLIFY = FALSE) %>%
  lapply(filter, log2FoldChange < 0)

names(targets50Top2) <- names(deseq)

saveRDS(targets50Top2, file = "/Users/effieklimi/Documents/novel-mirna/results/rds/vsmc-multimir.rds")
