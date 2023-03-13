library(tidyverse)
library(stringi)

setwd("/Users/effieklimi/Documents/novel-mirna/")

# GENCODE annotation:
gencodeV26 <- read.delim("gencode_v26_gtf_table.txt", header = FALSE, stringsAsFactors = FALSE)
annotation <- gencodeV26[, c(4, 7, 6, 1:3, 5)]
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")

vsmcRsem <- list.files(
  "rnaseq-data/rsem-vsmc",
  pattern = "*genes.results",
  full.names = TRUE
)[c(1,12,23,
2,13,24,
11,22,33,
10,21,32,
4,15,26,
5,16,27,
7,18,29,
9,20,31,
3,14,25,
6,17,28,
8,19,30)]

endosRsem <- list.files(
  "rnaseq-data/rsem-endos",
  pattern = "*genes.results",
  full.names = TRUE
  )[c(1,11,21,
2,12,22,
10,20,30,
4,14,24,
5,15,25,
7,17,27,
9,19,29,
3,13,23,
6,16,26,
8,18,28)]

vsmcNames <- gsub(
  "RSEM_JEI01_|-\\d[p]|.genes.results","",
        list.files(
          "rnaseq-data/rsem-vsmc",
          pattern = "*genes.results",
          full.names = FALSE
        )
  )

endosNames <- gsub(
  "RSEM_|-\\d[p]|_S\\d+|.genes.results","",
        list.files(
          "rnaseq-data/rsem-endos",
          pattern = "*genes.results",
          full.names = FALSE
        )
  )

vsmcENSEMBL <- read.table(vsmcRsem[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1] # access column 1 from file 1 - ENSEMBL IDs
# Count Table
vsmcCount <- do.call(cbind, lapply(vsmcRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 5]))
vsmcCountTable <- data.frame(vsmcENSEMBL, vsmcCount, stringsAsFactors = FALSE)
colnames(vsmcCountTable) <- c("ENSEMBL", vsmcNames)
vsmcCountTable <- vsmcCountTable[, c(1,2,13,24,3,14,25,12,23,34,11,22,33,5,16,27,6,17,28,8,19,30,10,21,32,4,15,26,7,18,29,9,20,31)]
vsmcCountTable <- merge(annotation, vsmcCountTable, by = 1)
# FPKM Table
vsmcFpkm <- do.call(cbind, lapply(vsmcRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 7]))
vsmcFpkmTable <- data.frame(vsmcENSEMBL, vsmcFpkm, stringsAsFactors = FALSE)
colnames(vsmcFpkmTable) <- c("ENSEMBL", vsmcNames)
vsmcFpkmTable <- vsmcFpkmTable[, c(1,2,13,24,3,14,25,12,23,34,11,22,33,5,16,27,6,17,28,8,19,30,10,21,32,4,15,26,7,18,29,9,20,31)]
vsmcFpkmTable <- merge(annotation, vsmcFpkmTable, by = 1)

write.csv(vsmcCountTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcCounts.csv", row.names = FALSE)
write.csv(vsmcFpkmTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcFpkm.csv", row.names = FALSE)


endosENSEMBL <- read.table(endosRsem[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1] # access column 1 from file 1 - ENSEMBL IDs
# Count Table
endosCount <- do.call(cbind, lapply(endosRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 5]))
endosCountTable <- data.frame(endosENSEMBL, endosCount, stringsAsFactors = FALSE)
colnames(endosCountTable) <- c("ENSEMBL", endosNames)
endosCountTable <- endosCountTable[, c(1,2,12,22,3,13,23,11,21,31,5,15,25,6,16,26,8,18,28,10,20,30,4,14,24,7,17,27,9,19,29 )]
endosCountTable <- merge(annotation, endosCountTable, by = 1)
# FPKM Table
endosFpkm <- do.call(cbind, lapply(endosRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 7]))
endosFpkmTable <- data.frame(endosENSEMBL, endosFpkm, stringsAsFactors = FALSE)
colnames(endosFpkmTable) <- c("ENSEMBL", endosNames)
endosFpkmTable <- endosFpkmTable[, c(1,2,12,22,3,13,23,11,21,31,5,15,25,6,16,26,8,18,28,10,20,30,4,14,24,7,17,27,9,19,29 )]
endosFpkmTable <- merge(annotation, endosFpkmTable, by = 1)

write.csv(endosCountTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endosCounts.csv", row.names = FALSE)
write.csv(endosFpkmTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endosFpkm.csv", row.names = FALSE)