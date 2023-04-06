library("tidyverse")
library("reshape2")


setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
expressed <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")

# get all DE genes 
vsmcDe <- 
  readRDS("results/rds/vsmc-deseq2.rds") %>% 
  lapply(as_tibble) %>%
  lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
  lapply(filter, name %in% expressed$name)

endosDe <- 
  readRDS("results/rds/endos-deseq2.rds") %>% 
  lapply(as_tibble) %>%
  lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
  lapply(filter, name %in% expressed$name)

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
multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)

names(multimir) <- miRNAnames
multimir <- multimir[c(2, 3, 5, 7, 1, 4, 6)]
targets <- do.call(rbind, multimir)
genes <-
    split(targets, duplicated(targets$name) | duplicated(targets$name, fromLast = TRUE)) %>%
    lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
    lapply(data.frame)

duplicates <- genes$`TRUE`$name
unique <- genes$`FALSE`$name

vsmc <- data.frame(
    multimirUnique = c(
        length(multimir[multimir[[2]]$name %in% unique]),
        length(multimir[multimir[[3]]$name %in% unique]),
        length(multimir[multimir[[5]]$name %in% unique]),
        length(multimir[multimir[[7]]$name %in% unique]),
        length(multimir[multimir[[1]]$name %in% unique]),
        length(multimir[multimir[[4]]$name %in% unique]),
        length(multimir[multimir[[6]]$name %in% unique])
    ),
    multimirDupl = c(
        length(multimir[multimir[[2]]$name %in% duplicates]),
        length(multimir[multimir[[3]]$name %in% duplicates]),
        length(multimir[multimir[[5]]$name %in% duplicates]),
        length(multimir[multimir[[7]]$name %in% duplicates]),
        length(multimir[multimir[[1]]$name %in% duplicates]),
        length(multimir[multimir[[4]]$name %in% duplicates]),
        length(multimir[multimir[[6]]$name %in% duplicates])
    ),
    DRgenes = c(
        + nrow(vsmcDr[[2]]) - nrow(multimir[[2]]),
        + nrow(vsmcDr[[3]]) - nrow(multimir[[3]]),
        + nrow(vsmcDr[[5]]) - nrow(multimir[[5]]),
        + nrow(vsmcDr[[7]]) - nrow(multimir[[7]]),
        + nrow(vsmcDr[[1]]) - nrow(multimir[[1]]),
        + nrow(vsmcDr[[4]]) - nrow(multimir[[4]]),
        + nrow(vsmcDr[[6]]) - nrow(multimir[[6]])
    )
  )


rownames(vsmc) <- c(
  "miR-323a-3p", 
  "miR-449b-5p", 
  "miR-491-3p", 
  "miR-892b", 
  "miR-1827", 
  "miR-4774-3p", 
  "miR-5681b"
)
  
vsmcMelt <- as.data.frame(vsmc) %>%
rownames_to_column("sample") %>%
mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
melt()

pdf(file = "results/figures/vsmc-multimir-genenumber.pdf", width = 10, height = 6)
ggplot(vsmcMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = .8) +
  scale_fill_manual(labels = c("Shared targets", "Unique targets", "Not-predicted"), values = c("#325791", "#7797d1", "#a64086")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of Genes") +
  labs(fill = " ") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"))
dev.off()


endos50Top2 <-
  readRDS("results/rds/endos-multimir.rds") %>%
  lapply("[", , c(3, 2)) %>%
  lapply(as_tibble)