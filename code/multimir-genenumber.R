library("tidyverse")


setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

# get all DE genes 
vsmcDe <- readRDS("results/rds/p01/vsmc-deseq2-p01.rds") %>% lapply(as_tibble)
endosDe <- readRDS("results/rds/p01/endos-deseq2-p01.rds") %>% lapply(as_tibble)

# get all DR genes
vsmcDr <-
  lapply(vsmcDe, "[", , c(7, 3)) %>%
  lapply(as_tibble) %>%
  lapply(filter, log2FoldChange < 0)

endosDr <-
  lapply(endosDe, "[", , c(7, 3)) %>%
  lapply(as_tibble) %>%
  lapply(filter, log2FoldChange < 0)

# get all UR genes
vsmcUr <-
  lapply(vsmcDe, "[", , c(7, 3)) %>%
  lapply(as_tibble) %>%
  lapply(filter, log2FoldChange > 0)

endosUr <-
  lapply(endosDe, "[", , c(7, 3)) %>%
  lapply(as_tibble) %>%
  lapply(filter, log2FoldChange > 0)



# Get all mltimir candidates
vsmc50Top2 <-
  readRDS("results/rds/p01/targets50top2-vsmc-p01.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)


endos50Top2 <-
  readRDS("results/rds/p01/targets50top2-endos-p01.rds") %>%
  lapply("[", , c(3, 2)) %>%
  lapply(as_tibble)



vsmc <- data.frame(
    URgenes = c(
        nrow(vsmcUr[[2]]),
        nrow(vsmcUr[[3]]),
        nrow(vsmcUr[[5]]),
        nrow(vsmcUr[[7]]),
        nrow(vsmcUr[[1]]),
        nrow(vsmcUr[[4]]),
        nrow(vsmcUr[[6]])

    ),
    DRgenes = c(
        nrow(vsmcDr[[2]]),
        nrow(vsmcDr[[3]]),
        nrow(vsmcDr[[5]]),
        nrow(vsmcDr[[7]]),
        nrow(vsmcDr[[1]]),
        nrow(vsmcDr[[4]]),
        nrow(vsmcDr[[6]])
    ),
    multimir = c(
        nrow(vsmc50Top2[[2]]),
        nrow(vsmc50Top2[[3]]),
        nrow(vsmc50Top2[[5]]),
        nrow(vsmc50Top2[[7]]),
        nrow(vsmc50Top2[[1]]),
        nrow(vsmc50Top2[[4]]),
        nrow(vsmc50Top2[[6]])
    )
  )



endos <- data.frame(
    URgenes = c(
        nrow(endosUr[[2]]),
        nrow(endosUr[[3]]),
        nrow(endosUr[[5]]),
        nrow(endosUr[[7]]),
        nrow(endosUr[[1]]),
        nrow(endosUr[[4]]),
        nrow(endosUr[[6]])

    ),
    DRgenes = c(
        nrow(endosDr[[2]]),
        nrow(endosDr[[3]]),
        nrow(endosDr[[5]]),
        nrow(endosDr[[7]]),
        nrow(endosDr[[1]]),
        nrow(endosDr[[4]]),
        nrow(endosDr[[6]])
    ),
    multimir = c(
        nrow(endos50Top2[[2]]),
        nrow(endos50Top2[[3]]),
        nrow(endos50Top2[[5]]),
        nrow(endos50Top2[[7]]),
        nrow(endos50Top2[[1]]),
        nrow(endos50Top2[[4]]),
        nrow(endos50Top2[[6]])
    )
  )

  rownames(vsmc) <- c("hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-491-3p", "hsa-miR-892b", "hsa-miR-1827", "hsa-miR-4774-3p", "hsa-miR-5681b")
  rownames(endos) <- c("hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-491-3p", "hsa-miR-892b", "hsa-miR-1827", "hsa-miR-4774-3p", "hsa-miR-5681b")

  write.csv(vsmc, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmc-genenumber-p01.csv")
  write.csv(endos, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endos-genenumber-p01.csv")