library("DESeq2")
library("dplyr")
library("tibble")
library("purrr")
library(foreach)



setwd("/Users/effieklimi/Documents/novel-mirna/")
gencodeV26 <- read.delim(
  "gencode_v26_gtf_table.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
annotation <- gencodeV26[, c(4, 7, 6, 1:3, 5)]
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")


type <- c(rep("paired-end", 12))
patient <- factor(rep(c("p1", "p2", "p3"), 4))
condition <- factor(c(
  rep("fbs02", 3),
  rep("ip", 3),
  rep("mock", 3),
  rep("mirctrl", 3)
  )
)

# Count data + meta data:
countTable <- read.csv(
  "results/tables/vsmcCounts.csv",
  header = TRUE
)

countTableDESeq <- sapply(countTable[, c(8:19)], as.integer)
row.names(countTableDESeq) <- countTable[, 1]
metadata <- data.frame(
  row.names = colnames(countTableDESeq),
  condition = condition,
  type = type,
  patient = patient
)

# Making DESeqDataSet object - adding design (correting for "patient"):
dds <- DESeqDataSetFromMatrix(
  countData = countTableDESeq,
  colData = metadata,
  design = ~ patient + condition
)

# Setting FBS02 as ref (can't use "contrast" due to running apeglm shrinkage):
dds$condition <- relevel(dds$condition, ref = "fbs02")
# Running DESeq2 with a Wald test:
dds <- DESeq(dds, test = "Wald")

resParams <- lapply(resultsNames(dds)[c(4:6)], list)
deseqResults <-
  foreach(contrast = resParams) %do% {
    res <- 
      results(
        dds,
        contrast = contrast,
        independentFiltering = TRUE,
        pAdjustMethod = "BH", # default
        alpha = 0.05
      )
  }

  shrinkResults <-
  foreach(contrast = resParams, res = deseqResults) %do% {

    shrink <-
    lfcShrink(
        dds,
        res = res,
        contrast = contrast,
        type = "ashr",
        svalue = FALSE,
        returnList = FALSE,
        format = "DataFrame"
      ) %>%
        data.frame() %>%
        rownames_to_column(var = "EnsID") %>%
        as_tibble() %>%
        filter(padj < 0.05) %>%
        merge(annotation, by = 1, all.x = FALSE)

  }

names(shrinkResults) <- resultsNames(dds)[c(4:6)]

saveRDS(shrinkResults, "results/rds/vsmc-deseq2-controls.rds")