library("clusterProfiler")
library("tidyverse")
library("purrr")
library("org.Hs.eg.db")
library("scales")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
targets <- readRDS("results/rds/vsmc-multimir.rds")

fpkm <- read.csv(
  "results/tables/vsmcExpressed.csv",
  header = TRUE
) %>% as_tibble()

# For the pathway analysis we took the uppermost quartile wrt logfoldchange
geneLists50Top2BP <-
  lapply(targets, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe) %>%
  lapply(filter, value < quantile(value, .25)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO,
      universe      = fpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1) %>%
  lapply(function(x) x@result)