library(tidyverse)

setwd("/Users/effieklimi/Documents/novel-mirna/")

# GENCODE annotation:
gencode_v26_gft_table <-
  read.csv(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz",
    header = FALSE,
    stringsAsFactors = FALSE
  )

gencode_v26_gft_table <- read.delim(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVEC RNA sequencing/gencode_v26_gtf_table.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
annotation <- gencode_v26_gft_table[, c(4, 7, 6, 1:3, 5)]
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")
#################################

######## read RSEM files and extract count and FPKM #########
#### CHANGE PATH to the RSEM folder you want to analyse ####
# "//cmvm.datastore.ed.ac.uk/cmvm/scs/users/username/..."

rsemVsmc <- list.files(
  "rnaseq-data/rsem-vsms",
  pattern = "*genes.results",
  full.names = TRUE
)

rsemEndos <- list.files("rnaseq-data/rsem-endos", pattern = "*genes.results", full.names = TRUE)


filenamesShort <- list.files("/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVEC RNA sequencing/RSEM HSVEC", pattern = "*genes.results", full.names = FALSE)
sampleNames <- gsub(
  "RSEM_", "",
  gsub(
    "-\\d[p]", "",
    gsub(
      "R-", "R",
      gsub(
        "_S\\d+.genes.results", "",
        gsub("0-2.genes.results", "0-2", filenamesShort)
      )
    )
  )
)
ENSEMBL <- read.table(filenames[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1] # access column 1 from file 1 - ENSEMBL IDs
# Count Table
readCount <- do.call(cbind, lapply(filenames, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 5]))
countTable <- data.frame(ENSEMBL, readCount, stringsAsFactors = FALSE)
colnames(countTable) <- c("ENSEMBL", sampleNames)
countTable <- merge(annotation, countTable, by = 1)
# FPKM Table
fpkm <- do.call(cbind, lapply(filenames, function(fn) read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 7]))
fpkmTable <- data.frame(ENSEMBL, fpkm, stringsAsFactors = FALSE)
colnames(fpkmTable) <- c("ENSEMBL", sampleNames)
fpkmTable <- merge(annotation, fpkmTable, by = 1)

write.csv(countTable, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVEC RNA sequencing/HSVEC_miRNAOE_counts.csv", row.names = FALSE)
write.csv(fpkmTable, file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVEC RNA sequencing/HSVEC_miRNAOE_FPKM.csv", row.names = FALSE)

#############################################################################################
#############################################################################################
#############################################################################################
