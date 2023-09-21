library("clusterProfiler")
library("tidyverse")
library("org.Hs.eg.db")
library(GOSemSim)
library(rrvgo)

options(warn = 1)
setwd("/Users/effieklimi/Documents/phd/novel-mirna/")

mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
targets <- readRDS("results/rds/vsmc-multimir-gobp-unclustered-filtered-deseq2.rds")

simMatrix <-
    lapply(targets, "[", , 1) %>%
    map(~calculateSimMatrix(.x, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel"))
saveRDS(simMatrix, "results/rds/semmantic-simmilarity/semmantic-sim-gobp-matrix.rds")


# load the matrix:
simMatrix <- readRDS("results/rds/semmantic-simmilarity/semmantic-sim-gobp-matrix.rds")

score <- lapply(targets, "[", , 5) %>% lapply(., function(item) -log10(item))

scoresLists <-
    score %>%
    map2(lapply(targets, "[", , 1), setNames)

reducedTerms <-
    simMatrix %>%
    map2(scoresLists, ~reduceSimMatrix(.x, .y, orgdb = "org.Hs.eg.db", threshold = 0.7))

pdf(file = "results/figures/semantic-similarity/semmanticsim-gobp-scatterplot.pdf", width = 8, height = 7)
map2(simMatrix, reducedTerms, scatterPlot)
dev.off()

pdf(file = "results/figures/semantic-similarity/semmanticsim-gobp-scatterplot.pdf", width = 8, height = 7)
treemapPlot(reducedTerms)
dev.off()
