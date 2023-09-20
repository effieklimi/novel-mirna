
options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

ensIdSplit <- function(x) { strsplit(x, split = ".", fixed = TRUE)[[1]][1] }
mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")



multimirRaw <-
    readRDS("results/rds/vsmc-multimir-frequency.rds") %>%
    lapply(as_tibble) %>%
    lapply(tibble::deframe) %>%
    lapply(sort, decreasing = TRUE) %>%
    lapply(tibble::enframe)
names(multimirRaw) <- mirNames
