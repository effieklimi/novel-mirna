library(tidyverse)
library(reshape2)

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

viability <- read.csv(
  "mirna-screen-data/r2-viability.csv",
  header = TRUE
) %>% as_tibble() %>% na.omit


toxic <- filter(viability, Z_score < -1.65)

toxicZscore <- toxic[, c("ID", "Z_score")] %>%
    arrange(Z_score) %>%
    mutate(ID = factor(ID, levels = ID)) %>% # This trick update the factor levels
    melt()

ggplot(toxicZscore, aes(x = ID, y = value, fill = variable)) +
    geom_col(stat = "identity", width = .8) +
    scale_fill_manual(
        values = rep(c("#bd4d94", "#4c71b5"), times = c(22, 1))) +
    scale_y_reverse() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

toxicNormalised <- toxic[, c("ID", "MCavg_Norm_Cell_Count")] %>%
    arrange(MCavg_Norm_Cell_Count) %>%
    mutate(ID = factor(ID, levels = ID)) %>% # This trick update the factor levels
    melt()

ggplot(toxicNormalised, aes(x = ID, y = value, fill = variable)) +
    geom_col(stat = "identity", width = .8) +
    scale_fill_manual(
        values = rep(c("#bd4d94", "#4c71b5"), times = c(22, 1))) +
    scale_y_reverse() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))