library(ggplot2)
library(tidyverse)
library(reshape2)


vsmc <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/vsmc-genenumbers-down-p01.csv", row.names = 1)
endos <- read.csv("/Users/effieklimi/Documents/novel-mirna/results/tables/endos-genenumbers-down-p01.csv", row.names = 1)

vmscMelt <- t(vsmc) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
  melt()

endosMelt <- t(endos) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample = factor(sample, levels = sample)) %>% # This trick update the factor levels
  melt()



vsmcPlot <- ggplot(vmscMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e67cbc", "#2d518a", "#72bcea")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of DE Genes")

endosPlot <- ggplot(endosMelt, aes(y = value, x = sample, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ea66b1", "#72bcea")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Number of DE Genes")