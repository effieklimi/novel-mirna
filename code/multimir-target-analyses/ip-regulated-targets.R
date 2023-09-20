library("ComplexHeatmap")
library("tidyverse")
library("org.Hs.eg.db")
library("enrichR")
library("clusterProfiler")
library(dendextend)
library(plyr)

#
options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")
# naming lists
miRnames <- c("hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-491-3p", "hsa-miR-892b", "hsa-miR-1827", "hsa-miR-4774-3p", "hsa-miR-5681b")
colNames <- c( "miR-Control", "miR-Control", "miR-Control",
 "miR-323a-3p", "miR-323a-3p", "miR-323a-3p",
 "miR-449b-5p", "miR-449b-5p", "miR-449b-5p",
"miR-491-3p", "miR-491-3p", "miR-491-3p",
"miR-892b", "miR-892b", "miR-892b",
"miR-1827", "miR-1827", "miR-1827",
"miR-4774-3p", "miR-4774-3p", "miR-4774-3p",
"miR-5681b", "miR-5681b", "miR-5681b")
#

controls <-
  readRDS("results/rds/vsmc-deseq2-controls.rds")[[1]][, c(7, 3)] %>% 
  as_tibble() %>%
  filter(log2FoldChange > log2(1.5)) %>%
  tibble::deframe() %>%
  sort(decreasing = FALSE) %>%
  tibble::enframe()

multimir <-
  readRDS("results/rds/vsmc-multimir.rds")[c(2, 3, 5, 7, 1, 4, 6)] %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe)
names(multimir) <- miRnames
targets <- do.call(rbind, multimir)

regulatedTargets <- multimir %>% lapply(filter, name %in% controls$name)

# Portion of targets regulated by IP:
## ~~~ miR-323a-3p
mir323Perc <- data.frame(
  "miRNA" = c("miR-323a-3p", "miR-323a-3p"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[1]]) / nrow(multimir[[1]])),
       1 - (nrow(regulatedTargets[[1]]) / nrow(multimir[[1]]))
    )
)

## ~~~ miR-449b-5p
mir449Perc <- data.frame(
  "miRNA" = c("miR-449b-5p", "miR-449b-5p"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[2]]) / nrow(multimir[[2]])),
       1 - (nrow(regulatedTargets[[2]]) / nrow(multimir[[2]]))
    )
)


## ~~~ miR-491-3p
mir491Perc <- data.frame(
  "miRNA" = c("miR-491-3p", "miR-491-3p"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[3]]) / nrow(multimir[[3]])),
       1 - (nrow(regulatedTargets[[3]]) / nrow(multimir[[3]]))
    )
)


## ~~~ miR-892b
mir892bPerc <- data.frame(
  "miRNA" = c("miR-892b", "miR-892b"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[4]]) / nrow(multimir[[4]])),
       1 - (nrow(regulatedTargets[[4]]) / nrow(multimir[[4]]))
    )
)


## ~~~ miR-1827
mir1827Perc <- data.frame(
  "miRNA" = c("miR-1827", "miR-1827"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[5]]) / nrow(multimir[[5]])),
       1 - (nrow(regulatedTargets[[5]]) / nrow(multimir[[5]]))
    )
)


## ~~~ miR-4774-3p
mir4774Perc <- data.frame(
  "miRNA" = c("miR-4774-3p", "miR-4774-3p"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[6]]) / nrow(multimir[[6]])),
       1 - (nrow(regulatedTargets[[6]]) / nrow(multimir[[6]]))
    )
)


## ~~~ miR-5681b
mir5681bPerc <- data.frame(
  "miRNA" = c("miR-5681b", "miR-5681b"),
  "Clusters" = c("Up-regulated by IP", "Not up-regulated by IP"),
  "Portion" = c(
      (nrow(regulatedTargets[[7]]) / nrow(multimir[[7]])),
       1 - (nrow(regulatedTargets[[7]]) / nrow(multimir[[7]]))
    )
)

percentageIP <- rbind(mir323Perc, mir449Perc, mir491Perc, mir892bPerc, mir1827Perc, mir4774Perc, mir5681bPerc)
pdf(file = "results/figures/multimir-ip-regulation-percetages.pdf", width = 9, height = 4)
ggplot(percentageIP, aes(fill=Clusters, y=Portion, x = factor(miRNA, level=c('miR-323a-3p', 'miR-449b-5p', 'miR-491-3p', 'miR-892b', 'miR-1827', 'miR-4774-3p', 'miR-5681b')))) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = c("#d4d4d4", "#af6189")) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) +
    labs(fill = "Genes that are:") +
    xlab(' ') +
    ylab("Percentage of targets in each category")
dev.off()