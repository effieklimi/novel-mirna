library("DESeq2")
library("dplyr")
library("tibble")
library("purrr")


setwd("/Users/effieklimi/Documents/novel-mirna/")
gencodeV26 <- read.delim(
  "gencode_v26_gtf_table.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
annotation <- gencodeV26[, c(4, 7, 6, 1:3, 5)]
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")


type <- c(rep("paired-end", 30))
patient <- factor(rep(c("p1", "p2", "p3"), 10))
condition <- factor(c(
  rep("FBS02perc", 3),
  rep("FBS10perc", 3),
  rep("mirctrl", 3),
  rep("mir323", 3),
  rep("mir449b", 3),
  rep("mir491", 3),
  rep("mir892b", 3),
  rep("mir1827", 3),
  rep("mir4774", 3),
  rep("mir5681b", 3)
  )
)

# Count data + meta data:
countTable <- read.csv(
  "results/tables/endosCounts.csv",
  header = TRUE
)

countTableDESeq <- sapply(countTable[, c(8:37)], as.integer)
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
dds$condition <- relevel(dds$condition, ref = "mirctrl")
# Running DESeq2 with a Wald test:
dds <- DESeq(dds, test = "Wald")

resParams <- lapply(resultsNames(dds)[c(6:12)], list)
deseqResults <-
  foreach(contrast = resParams) %do% {

    res <- 
      results(
        dds,
        contrast = contrast,
        independentFiltering = TRUE,
        pAdjustMethod = "BH", # default
        alpha = 0.01
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
        filter(padj < 0.01) %>%
        merge(annotation, by = 1, all.x = FALSE)

  }

names(shrinkResults) <- resultsNames(dds)[c(6:12)]

saveRDS(shrinkResults, "results/tables/endos-deseq2-results.rds")