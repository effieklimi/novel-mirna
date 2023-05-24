library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)
theme_set(theme_pubr())

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

viability <- read.csv(
    "mirna-screen-data/TriesteScreen_miRCTRLFC_Viability.csv",
    header = TRUE
  ) %>% as_tibble() %>% na.omit


viability2 <- viability
# Hide all of the text labels.
viability2$ID <- ""
# Let's just label these items.
ix_label <- which(viability2$Categories == "Interesting")
viability2$ID[ix_label] <- viability$ID[ix_label]


ggplot(viability2, aes(x=MCavg_Norm_Proliferative, y=MCavg_Norm_Cell_Count, label = ID)) + 
  geom_point(alpha=.4, size=2, aes(color=factor(Categories))) +
  scale_color_manual(values = c('darkred', "darkgoldenrod3", '#999999', "deepskyblue4", "black")) +
  geom_vline(xintercept = 1, linetype="dashed", color = "black", size = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0.637, linetype = "dashed", color = "black", size = 0.2, alpha = 0.7) +
  geom_label_repel(
    max.overlaps = Inf,
    segment.linetype = 2) +
  labs(title = "High-throughput screen results", subtitle = "Viability & fold change in proliferation per miRNA", caption = "(2046 miRNAs)") +
  xlab("Fold change (VS miR-CTRL mean)") +
  ylab("Viability (z-score of cell count)") +
  theme(legend.position = "bottom") +
  theme_classic()

Screen_viability <- melt(Screen_viability)
ggplot(Screen_viability, aes( y=value, x=variable, fill=variable)) + 
  see::geom_violinhalf(alpha=0.3, fill="#F398C5", color="#CF3B85") +
  labs(x = NULL) +
  scale_y_continuous(breaks=seq(0, 6, 1))
theme_pubclean()

screen_table <- read.csv("/Users/effieklimi/Documents/PhD/Novel miRs paper/TriesteScreen/TriesteScreen_FC_MC.csv", header = TRUE, sep = ",")
MCmean <- screen_table[,c(1,2)]
MCmean <- melt(MCmean)

MC1 <- screen_table[,c(1,3)]
MC1 <- melt(MC1)

MC2 <- screen_table[,c(1,4)]
MC2 <- melt(MC2)

MC3 <- screen_table[,c(1,5)]
MC3 <- melt(MC3)

MC4 <- screen_table[,c(1,6)]
MC4 <- melt(MC4)

screen_table <- melt(screen_table)



ggplot(MCmean, aes(log2(value), fill = variable, colour = variable)) +
  geom_density(alpha=0.1)



ggplot(MCmean, aes( y=value, x=variable, fill=variable)) + 
  see::geom_violinhalf(alpha=0.3, position = position_dodge(width = 3),size=1,color=NA, fill="lightcoral", width=0.7) +
  geom_point(size=0.5, position = position_jitterdodge(), alpha=0.2) + 
  ggbeeswarm::geom_quasirandom(shape = 20,size=0.5, dodge.width = .75, color="blue",alpha=.5,show.legend = F)+
  scale_color_manual(values=c('lightcoral ', 'azure3', "lightcoral","black" )) +
  geom_boxplot(notch = FALSE, width=0.1, outlier.size = -1, color="black",lwd=0.5, alpha = 0.7, fill="lightcoral") +
  labs(x = NULL) +
  theme_minimal()+
  ylab(  c("Fold proliferation")  )  +
  xlab(  c(" ")  )  +
  rremove("legend.title")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.line = element_line(colour = "black",size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = c(0.95, 0.85))+
  font("xylab",size=15)+  
  font("xy",size=15)+ 
  font("xy.text", size = 15) +  
  font("legend.text",size = 15)+
  guides(fill = guide_legend(override.aes = list(alpha = 1,color="black")))
  theme_pubclean()



ggplot(MCmean, aes( y=value, x=variable, fill=variable)) + 
  see::geom_violinhalf(alpha=0.3, fill="#F398C5", color="#CF3B85") +
  labs(x = NULL) + 
  coord_flip() +
  scale_y_continuous(breaks=seq(0,6,1))
  theme_pubclean()
  
# geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.2,alpha=0.5) +
  