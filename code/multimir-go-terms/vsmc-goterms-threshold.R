library("clusterProfiler")
library("tidyverse")
library("purrr")
library("org.Hs.eg.db")
library("scales")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

mirNames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
targets <- readRDS("results/rds/vsmc-multimir-threshold.rds")

fpkm <- read.csv(
  "results/tables/vsmcExpressed.csv",
  header = TRUE
) %>% as_tibble()

geneLists50Top2BP <-
  lapply(targets, "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(names) %>%
  map(enrichGO,
      universe      = fpkm$name,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff = 1,
      minGSSize = 10)

bpEnrichTable <- lapply(geneLists50Top2BP, function(x) x@result)
saveRDS(geneLists50Top2BP, file = "results/rds/vsmc-multimir-gobp-unclustered-filtered-deseq2.rds")

bpSimplify <- lapply(bpEnrichment, simplify, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
bpSimpTable <- lapply(bpSimplify, function(x) x@result)

write.table(geneLists50Top2BP[[1]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir1827.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[2]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir323.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[3]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir449b.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[4]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir4774.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[5]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir491.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[6]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir5681b.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
write.table(geneLists50Top2BP[[7]][1:50, c(1, 5)], file = "/Users/effieklimi/Documents/novel-mirna/results/tables/goterm-semantic-similarity/unclustered-filtered-deseq2/mir892b.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("% GOterm", "enrichment_P-value"))
