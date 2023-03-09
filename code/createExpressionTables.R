library(tidyverse)

setwd("/Users/effieklimi/Documents/novel-mirna/")

# GENCODE annotation:
gencodeV26 <- read.delim("gencode_v26_gtf_table.txt", header = FALSE, stringsAsFactors = FALSE)
annotation <- data.frame(
  ensembl = gencodeV26$V4,
  name = gencodeV26$V7,
  type = gencodeV26$V6,
  chr = gencodeV26$V1,
  start = gencodeV26$V2,
  end = gencodeV26$V3,
  str = gencodeV26$V5
)
# can also use the following two lines: annotation <- gencodeAnnot[ , c(4, 7, 6, 1:3, 5)]
#to do the same as with the data.frame function above: colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")
#################################

######## read RSEM files and extract count and FPKM #########
#### CHANGE PATH to the RSEM folder you want to analyse ####
# "//cmvm.datastore.ed.ac.uk/cmvm/scs/users/username/..."

rsemVsmc <- list.files(
  "rnaseq-data/rsem-vsmc",
  pattern = "*genes.results",
  full.names = TRUE
)

rsemEndos <- list.files(
  "rnaseq-data/rsem-endos",
  pattern = "*genes.results",
  full.names = TRUE
  )

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
