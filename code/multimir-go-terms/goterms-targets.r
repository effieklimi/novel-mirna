library("org.Hs.eg.db")
library("tidyverse")


deseq <- readRDS("results/rds/vsmc-deseq2.rds")
genes <- 
  lapply(deseq, "[", , c(1, 3, 7)) %>%
  map(~ mutate(.x, EnsID = ensIdSplit(EnsID))) %>%
  lapply(as_tibble) %>%
  lapply(filter, name %in% fpkm$name)
names(genes) <- mirNames


multimir <-
  readRDS("results/rds/vsmc-multimir.rds") %>%
  lapply( "[", , c(3, 2)) %>%
  lapply("as_tibble") %>%
  lapply(tibble::deframe) %>%
  lapply(sort, decreasing = FALSE) %>%
  lapply(tibble::enframe)


cellcycle <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", columns="SYMBOL")
cellcycleTargets <- multimir %>%
    lapply(filter, name %in% cellcycle$SYMBOL)

motility <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0016477", columns="SYMBOL")
motilityTargets <- multimir %>%
    lapply(filter, name %in% motility$SYMBOL)

cationtrans <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0006812", columns="SYMBOL")
cationtransTargets <- multimir %>%
    lapply(filter, name %in% cationtrans$SYMBOL)

proteinmod <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0006486", columns="SYMBOL")
proteinmodTargets <- multimir %>%
    lapply(filter, name %in% proteinmod$SYMBOL)

lipidtrans <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0006869", columns="SYMBOL")
lipidtransTargets <- multimir %>%
    lapply(filter, name %in% lipidtrans$SYMBOL)

gtpase <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007264", columns="SYMBOL")
gtpaseTargets <- multimir %>%
    lapply(filter, name %in% gtpase$SYMBOL)

vascdev <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0001944", columns="SYMBOL")
vascdevTargets <- multimir %>%
    lapply(filter, name %in% vascdev$SYMBOL)

leukochem <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0030595", columns="SYMBOL")
leukochemTargets <- multimir %>%
    lapply(filter, name %in% leukochem$SYMBOL)

inflam <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0006954", columns="SYMBOL")
inflamTargets <- multimir %>%
    lapply(filter, name %in% inflam$SYMBOL)


    





# cell cycle graphs - miR-323a, miR-1827, miR-449b, miR-491, miR-892b
pdf("results/figures/cellcycletargets/cellcycle-miR1827.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[1]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-1827") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR323.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[2]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-323a-3p") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR449.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[3]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-449b-5p") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR4774.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[4]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-4774-3p") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR491.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[5]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-491-3p") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR5681.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[6]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-5681b") +
  theme_1()
dev.off()

pdf("results/figures/cellcycletargets/cellcycle-miR892b.pdf", width = 4.5, height = 4)
ggplot(cellcycleTargets[[7]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#7192ce") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-892b") +
  theme_1()
dev.off()


# motility graphs - miR-4774-3p, miR-491
pdf("results/figures/migrationtargets/migration-miR4774.pdf", width = 4.5, height = 4)
ggplot(motilityTargets[[4]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#e28e8e") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-4774-3p") +
  theme_1()
dev.off()

pdf("results/figures/migrationtargets/migration-miR491.pdf", width = 4.5, height = 4)
ggplot(motilityTargets[[5]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#e28e8e") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-491-3p") +
  theme_1()
dev.off()

# glyco graphs - miR-4774-3p, miR-5681b
pdf("results/figures/proteinmodTargets/proteinmod-miR4774.pdf", width = 4.5, height = 4)
ggplot(proteinmodTargets[[4]], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#e2ce8e") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-4774-3p") +
  theme_1()
dev.off()

pdf("results/figures/proteinmodTargets/proteinmod-miR5681.pdf", width = 4.5, height = 4)
ggplot(proteinmodTargets[[6]], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#e2ce8e") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-5681b") +
  theme_1()
dev.off()

# ion transport (calcium potassium) graphs - miR-449b miR-491
pdf("results/figures/iontransTargets/iontrans-miR449.pdf", width = 4.5, height = 4)
ggplot(iontransTargets[[3]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#c87ebf") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-449b-5p") +
  theme_1()
dev.off()

pdf("results/figures/iontransTargets/iontrans-miR491.pdf", width = 4.5, height = 4)
ggplot(iontransTargets[[5]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#c87ebf") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-491-3p") +
  theme_1()
dev.off()

# vascular development graphs - miR-892b
pdf("results/figures/vascdevTargets/vascdev-miR892b.pdf", width = 4.5, height = 4)
ggplot(vascdevTargets[[7]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#aeaeae") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-892b") +
  theme_1()
dev.off()

# leuko
pdf("results/figures/leukochemTargets/leukochemTargets-miR5681b.pdf", width = 4.5, height = 4)
ggplot(leukochemTargets[[6]], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#aeaeae") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-5681b") +
  theme_1()
dev.off()

# lipid 
pdf("results/figures/lipidtransTargets/lipidtransTargets-miR491.pdf", width = 4.5, height = 4)
ggplot(lipidtransTargets[[5]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#9ec87e") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-491-3p") +
  theme_1()
dev.off()

# inflamTargets
pdf("results/figures/inflamTargets/inflam-miR491.pdf", width = 4.5, height = 4)
ggplot(inflamTargets[[5]][1:10,], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#ab8ee2") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-491-3p") +
  theme_1()
dev.off()

pdf("results/figures/inflamTargets/inflam-miR5681.pdf", width = 4.5, height = 4)
ggplot(inflamTargets[[6]], aes(x = value, y = reorder(name, -value))) +
  geom_col(fill = "#ab8ee2") +
  scale_x_reverse() +
  labs(x = "Log Fold Change", y = " ", title = "miR-5681b") +
  theme_1()
dev.off()




# gtpase




























theme_1 <- function(){
    theme_bw() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -3),
        panel.background = element_rect(fill = "transparent", colour = NA),
        #panel.border = element_blank(),
        panel.grid.major = element_line(colour="grey100", size=0.1),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
        plot.title = element_text(size = 15, vjust = 1, hjust = 0),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_rect(
            color = "black",
            fill = "transparent",
            size = 3, linetype = "blank")
    ) 
}