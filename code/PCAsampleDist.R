library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(svglite)
library(DESeq2)


# PCA #
countTable <- read.csv("/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/HSVSMC_miRNAOE_counts.csv",
  header = TRUE, stringsAsFactors = FALSE
)
#### Group columns by condition (arithmetically)
countTablePCA <- sapply(
  countTable[, c(
    8, 19, 30,  # 02
    9, 20, 31,  # IP
    18, 29, 40, # Mock
    17, 28, 39, # miRCTRL
    11, 22, 33, # miR323
    12, 23, 34, # miR449b
    14, 25, 36, # miR491
    16, 27, 38, # miR892b
    10, 21, 32, # miR1827
    13, 24, 35, # miR4774
    15, 26, 37  # miR5681b
  )], 
  as.integer
)
row.names(countTable) <- countTable[, 1]
# modify conditions and number of samples:
type <- c(rep("paired-end", 33))
patient <- rep(c("p1", "p2", "p3"), 11)
condition <- c(
  rep("Quiescent", 3),
  rep("IL-1A/PDGF", 3),
  rep("Mock", 3),
  rep("miR-Control", 3),
  rep("miR-323a-3p", 3),
  rep("miR-449b-5p", 3),
  rep("miR-491-3p", 3),
  rep("miR-892b", 3),
  rep("miR-1827", 3),
  rep("miR-4774-3p", 3),
  rep("miR-5681b", 3)
)
PCAannotCountData <- data.frame(condition = condition, type = type, patient = patient)
row.names(PCAannotCountData) <- colnames(countTablePCA)
# DESeq2: preparing matrix and performing analysis
ddsPCA <- DESeqDataSetFromMatrix(
  countData = countTablePCA,
  colData = PCAannotCountData,
  design = ~ patient + condition
)
ddsPCA <- DESeq(ddsPCA)
rldPCA <- rlog(ddsPCA, blind = FALSE)
# Batch correction with limma package
assay(rldPCA) <- limma::removeBatchEffect(assay(rldPCA),
  batch = rldPCA$patient,
  design = model.matrix(~condition, colData(rldPCA))
)
# PCA plot
PCA <- plotPCA(rldPCA, intgroup = c("condition"), ntop = 500, returnData = TRUE)
rldPCAnames <- rownames(colData(rldPCA))
percentVar <- round(100 * attr(PCA, "percentVar"))
PCA$Patient <- patient
PCA$Condition <- condition
colnames(PCA) <- c("PC1", "PC2", "group", "Condition", "name", "Patient")

#### theme #####
theme_1 <- function() {
  theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.title.y = element_text(vjust = +3),
      axis.title.x = element_text(vjust = -3),
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_line(colour = "grey100", linewidth = 0.1),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 15, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.key = element_blank(),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        linewidth = 3, linetype = "blank"
      )
    )
}
# Theme: #
################

# PCA #################################################################################################################################
PCAplot <- ggplot(PCA, aes(PC1, PC2, color = Condition, shape = Patient, frame = TRUE)) + #
  geom_point(size = 5) + #
  xlab(paste0("Principal Component 1: ", percentVar[1], "% variance")) + #
  ylab(paste0("Principal Component 2: ", percentVar[2], "% variance")) + #
  # geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  scale_color_manual(
    values = c(
      "grey", # Quiescent condition
      "#B7154B", "#E52464", "#F07FA5", # IP-stimulated controls
      "#0D1E30", "#1E4671", "#2B64A1", "#3E82CC", "#5E97D4", "#8EB6E1", "#CFE0F2" # miRNA OE
    ), 
    breaks = c(
      "Quiescent",
      "IL-1A/PDGF", "Mock", "miR-Control",
      "miR-323a-3p", "miR-449b-5p", "miR-491-3p", "miR-892b", "miR-1827", "miR-4774-3p", "miR-5681b"
    )
  ) +
  scale_x_continuous(limits = c(-40, 17)) +
  theme_1()
#######################################################################################################################################


# Sample distances #################################################
sampleDists <- dist(t(assay(rldPCA))) #
sampleDistMatrix <- as.matrix(sampleDists) #
rownames(sampleDistMatrix) <- paste(rldPCA$condition) #
colnames(sampleDistMatrix) <- paste(rldPCA$condition) #
#
colors <- colorRampPalette(rev(brewer.pal(9, "BuPu")))(255) #
#
distPheatmap <- pheatmap(sampleDistMatrix, #
  clustering_distance_rows = sampleDists, #
  clustering_distance_cols = sampleDists, #
  col                      = colors, #
  fontsize_row             = 13, #
  fontsize_col             = 13
) #
####################################################################


# Save figures
setwd("/Users/effieklimi/Documents/PhD/miRNA screening paper/Figures/")
ggsave(
  filename = "sampleDistances.svg",
  plot = distPheatmap,
  width = 2500,
  height = 2500,
  dpi = 300,
  units = "px",
  device = "svg"
)

ggsave(
  filename = "PCAplot.svg",
  plot = PCAplot,
  width = 2500,
  height = 2000,
  dpi = 300,
  units = "px",
  device = "svg"
)
