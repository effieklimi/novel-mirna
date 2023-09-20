library("tidyverse")
library("reshape2")


setwd("/Users/effieklimi/Documents/novel-mirna/")
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
expressed <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmcExpressed.csv")

# get all DE genes 
vsmcDe <- 
  readRDS("results/rds/vsmc-deseq2-mimics.rds") %>% 
  lapply(as_tibble) %>%
  lapply(distinct, name, .keep_all = TRUE) %>% # remove duplicates
  lapply(filter, name %in% expressed$name)

endosDe <- 
  readRDS("results/rds/endos-deseq2-mimics.rds") %>% 
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
vsmc50Top2 <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply(as_tibble)


endos50Top2 <-
  readRDS("results/rds/endos-multimir.rds") %>%
  lapply("[", , c(3, 2)) %>%
  lapply(as_tibble)



vsmc <- data.frame(
    URgenes = c(
        nrow(vsmcUr[[2]]),
        nrow(vsmcUr[[3]]),
        nrow(vsmcUr[[5]]),
        nrow(vsmcUr[[7]]),
        nrow(vsmcUr[[1]]),
        nrow(vsmcUr[[4]]),
        nrow(vsmcUr[[6]])

    ),
    multimir = c(
        - nrow(vsmc50Top2[[2]]),
        - nrow(vsmc50Top2[[3]]),
        - nrow(vsmc50Top2[[5]]),
        - nrow(vsmc50Top2[[7]]),
        - nrow(vsmc50Top2[[1]]),
        - nrow(vsmc50Top2[[4]]),
        - nrow(vsmc50Top2[[6]])
    ),
    DRgenes = c(
        - nrow(vsmcDr[[2]]) + nrow(vsmc50Top2[[2]]),
        - nrow(vsmcDr[[3]]) + nrow(vsmc50Top2[[3]]),
        - nrow(vsmcDr[[5]]) + nrow(vsmc50Top2[[5]]),
        - nrow(vsmcDr[[7]]) + nrow(vsmc50Top2[[7]]),
        - nrow(vsmcDr[[1]]) + nrow(vsmc50Top2[[1]]),
        - nrow(vsmcDr[[4]]) + nrow(vsmc50Top2[[4]]),
        - nrow(vsmcDr[[6]]) + nrow(vsmc50Top2[[6]])
    )
  )



endos <- data.frame(
    URgenes = c(
        nrow(endosUr[[2]]),
        nrow(endosUr[[3]]),
        nrow(endosUr[[5]]),
        nrow(endosUr[[7]]),
        nrow(endosUr[[1]]),
        nrow(endosUr[[4]]),
        nrow(endosUr[[6]])

    ),
    multimir = c(
        -nrow(endos50Top2[[2]]),
        -nrow(endos50Top2[[3]]),
        -nrow(endos50Top2[[5]]),
        -nrow(endos50Top2[[7]]),
        -nrow(endos50Top2[[1]]),
        -nrow(endos50Top2[[4]]),
        -nrow(endos50Top2[[6]])
    ),
    DRgenes = c(
        - nrow(endosDr[[2]]) + nrow(endos50Top2[[2]]),
        - nrow(endosDr[[3]]) + nrow(endos50Top2[[3]]),
        - nrow(endosDr[[5]]) + nrow(endos50Top2[[5]]),
        - nrow(endosDr[[7]]) + nrow(endos50Top2[[7]]),
        - nrow(endosDr[[1]]) + nrow(endos50Top2[[1]]),
        - nrow(endosDr[[4]]) + nrow(endos50Top2[[4]]),
        - nrow(endosDr[[6]]) + nrow(endos50Top2[[6]])
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

pdf(file = "results/figures/vsmc-down-multimir-bar-complete.pdf", width = 8, height = 4.5)
ggplot(vsmcMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = .8) +
  scale_fill_manual(labels = c("Up-regulated Genes", "Down-regulated, Predicted as targets", "Down-regulated, Not predicted as targets"), 
                    values = c( "#b03280","#325791", "#7797d1")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of Genes") +
  labs(fill = " ") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 15, face = "bold"))
dev.off()

  rownames(endos) <- c(
    "miR-323a-3p", 
    "miR-449b-5p", 
    "miR-491-3p", 
    "miR-892b", 
    "miR-1827", 
    "miR-4774-3p", 
    "miR-5681b"
  )

endosMelt <- as.data.frame(endos) %>%
rownames_to_column("sample") %>%
mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
melt()

pdf(file = "results/figures/endos-down-multimir-bar-complete.pdf", width = 10, height = 4.5)
ggplot(endosMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = .8) +
  scale_fill_manual(labels = c("Up-regulated Genes", "Down-regulated, Predicted as targets", "Down-regulated, Not predicted as targets"), values = c("#c63886", "#325791", "#7797d1")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of Genes") +
  labs(fill = " ") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 15, face = "bold"))
dev.off()

  rownames(vsmc) <- c("miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b")
  rownames(endos) <- c("miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b")

  write.csv(vsmc, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmc-genenumbers-comprehensive-p01.csv")
  write.csv(endos, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endos-genenumbers-comprehensive-p01.csv")
















vsmc <- data.frame(
    URgenes = c(
        nrow(vsmcUr[[2]]),
        nrow(vsmcUr[[3]]),
        nrow(vsmcUr[[5]]),
        nrow(vsmcUr[[7]]),
        nrow(vsmcUr[[1]]),
        nrow(vsmcUr[[4]]),
        nrow(vsmcUr[[6]])

    ),
    DRgenes = c(
        nrow(vsmcDr[[2]]),
        nrow(vsmcDr[[3]]),
        nrow(vsmcDr[[5]]),
        nrow(vsmcDr[[7]]),
        nrow(vsmcDr[[1]]),
        nrow(vsmcDr[[4]]),
        nrow(vsmcDr[[6]])
    )
  )



endos <- data.frame(
    URgenes = c(
        nrow(endosUr[[2]]),
        nrow(endosUr[[3]]),
        nrow(endosUr[[5]]),
        nrow(endosUr[[7]]),
        nrow(endosUr[[1]]),
        nrow(endosUr[[4]]),
        nrow(endosUr[[6]])

    ),
    DRgenes = c(
        nrow(endosDr[[2]]),
        nrow(endosDr[[3]]),
        nrow(endosDr[[5]]),
        nrow(endosDr[[7]]),
        nrow(endosDr[[1]]),
        nrow(endosDr[[4]]),
        nrow(endosDr[[6]])
    )
  )

rownames(vsmc) <- c("miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b")
rownames(endos) <- c("miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b")

write.csv(vsmc, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmc-genenumbers-updown-p01.csv")
write.csv(endos, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endos-genenumbers-updown-p01.csv")




vsmc <- data.frame(
    multimir = c(
        nrow(vsmc50Top2[[2]]),
        nrow(vsmc50Top2[[3]]),
        nrow(vsmc50Top2[[5]]),
        nrow(vsmc50Top2[[7]]),
        nrow(vsmc50Top2[[1]]),
        nrow(vsmc50Top2[[4]]),
        nrow(vsmc50Top2[[6]])
    ),
    DRgenes = c(
        + nrow(vsmcDr[[2]]) - nrow(vsmc50Top2[[2]]),
        + nrow(vsmcDr[[3]]) - nrow(vsmc50Top2[[3]]),
        + nrow(vsmcDr[[5]]) - nrow(vsmc50Top2[[5]]),
        + nrow(vsmcDr[[7]]) - nrow(vsmc50Top2[[7]]),
        + nrow(vsmcDr[[1]]) - nrow(vsmc50Top2[[1]]),
        + nrow(vsmcDr[[4]]) - nrow(vsmc50Top2[[4]]),
        + nrow(vsmcDr[[6]]) - nrow(vsmc50Top2[[6]])
    )
  )


endos <- data.frame(
  multimir = c(
    nrow(endos50Top2[[2]]),
    nrow(endos50Top2[[3]]),
    nrow(endos50Top2[[5]]),
    nrow(endos50Top2[[7]]),
    nrow(endos50Top2[[1]]),
    nrow(endos50Top2[[4]]),
    nrow(endos50Top2[[6]])
  ),
  DRgenes = c(
    + nrow(endosDr[[2]]) - nrow(endos50Top2[[2]]),
    + nrow(endosDr[[3]]) - nrow(endos50Top2[[3]]),
    + nrow(endosDr[[5]]) - nrow(endos50Top2[[5]]),
    + nrow(endosDr[[7]]) - nrow(endos50Top2[[7]]),
    + nrow(endosDr[[1]]) - nrow(endos50Top2[[1]]),
    + nrow(endosDr[[4]]) - nrow(endos50Top2[[4]]),
    + nrow(endosDr[[6]]) - nrow(endos50Top2[[6]])
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

pdf(file = "results/figures/vsmc-down-multimir-bar.pdf", width = 8, height = 4.5)
ggplot(vsmcMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = .8) +
  scale_fill_manual(labels = c("Predicted as targets", "Non-predicted as targets"), values = c("#325791", "#7797d1")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of Genes") +
  labs(fill = " ") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 15, face = "bold"))
dev.off()

  rownames(endos) <- c(
    "miR-323a-3p", 
    "miR-449b-5p", 
    "miR-491-3p", 
    "miR-892b", 
    "miR-1827", 
    "miR-4774-3p", 
    "miR-5681b"
  )

endosMelt <- as.data.frame(endos) %>%
rownames_to_column("sample") %>%
mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
melt()

pdf(file = "results/figures/endos-down-multimir-bar.pdf", width = 8, height = 4.5)
ggplot(endosMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = .8) +
  scale_fill_manual(labels = c("Predicted as targets", "Non-predicted as targets"), values = c("#325791", "#7797d1")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of Genes") +
  labs(fill = " ") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 15, face = "bold"))
dev.off()

write.csv(vsmc, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/vsmc-genenumbers-down-p01.csv")
write.csv(endos, file = "/Users/effieklimi/Documents/novel-mirna/results/tables/endos-genenumbers-down-p01.csv")


