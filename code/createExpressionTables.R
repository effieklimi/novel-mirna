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
)

endosRsem <- list.files(
  "rnaseq-data/rsem-endos",
  pattern = "*genes.results",
  full.names = TRUE
  )

endosNames <- gsub(
  "RSEM_|-\\d[p]|_S\\d+|.genes.results","",
        list.files(
        "rnaseq-data/rsem-endos",
        pattern = "*genes.results",
        full.names = FALSE
        )
  )

vsmcNames <- gsub(
  "RSEM_JEI01_|-\\d[p]|.genes.results","",
        list.files(
        "rnaseq-data/rsem-vsmc",
        pattern = "*genes.results",
        full.names = FALSE
        )
  )


vsmcENSEMBL <- read.table(vsmcRsem[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1] # access column 1 from file 1 - ENSEMBL IDs
# Count Table
vsmcCount <- do.call(cbind, lapply(vsmcRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 5]))
vsmcCountTable <- data.frame(ENSEMBL, readCount, stringsAsFactors = FALSE)
colnames(vsmcCountTable) <- c("ENSEMBL", vsmcNames)
vsmcCountTable <- merge(annotation, vsmcCountTable, by = 1)
# FPKM Table
vsmcFpkm <- do.call(cbind, lapply(vsmcRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 7]))
vsmcFpkmTable <- data.frame(ENSEMBL, fpkm, stringsAsFactors = FALSE)
colnames(vsmcFpkmTable) <- c("ENSEMBL", vsmcNames)
vsmcFpkmTable <- merge(annotation, vsmcFpkmTable, by = 1)

write.csv(vsmcCountTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcCounts.csv", row.names = FALSE)
write.csv(vsmcFpkmTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcFpkm.csv", row.names = FALSE)


endosENSEMBL <- read.table(vsmcRsem[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1] # access column 1 from file 1 - ENSEMBL IDs
# Count Table
endosCount <- do.call(cbind, lapply(endosRsem, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 5]))
endosCountTable <- data.frame(ENSEMBL, endosCount, stringsAsFactors = FALSE)
colnames(endosCountTable) <- c("ENSEMBL", endosNames)
endosCountTable <- merge(annotation, endosCountTable, by = 1)
# FPKM Table
endosFpkmTable <- data.frame(ENSEMBL, endosFpkm, stringsAsFactors = FALSE)
colnames(endosFpkmTable) <- c("ENSEMBL", endosNames)
endosFpkmTable <- merge(annotation, endosFpkmTable, by = 1)

write.csv(endosCountTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endosCounts.csv", row.names = FALSE)
write.csv(endosFpkmTable, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endosFpkm.csv", row.names = FALSE)