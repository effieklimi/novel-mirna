library("DESeq2")
library("tidyverse")
library("tibble")
library("purrr")
library("foreach")
library("RColorBrewer")
library("pheatmap")

setwd("/Users/effieklimi/Documents/novel-mirna/")

gencodeV26 <- read.delim(
  "gencode_v26_gtf_table.txt",
  header = FALSE,
  stringsAsFactors = FALSE
);

annotation <- gencodeV26[, c(4, 7, 6, 1:3, 5)];
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str");


type <- c(rep("paired-end", 33));
patient <- factor(rep(c("p1", "p2", "p3"), 11));
condition <- factor(c(
    rep("Quiesced", 3),
    rep("IL-1A/PDGF-BB", 3),
    rep("Mock", 3),
    rep("miR-CTRL", 3),
    rep("miR-323", 3),
    rep("miR-449b", 3),
    rep("miR-491", 3),
    rep("miR-892b", 3),
    rep("miR-1827", 3),
    rep("miR-4774", 3),
    rep("miR-5681b", 3)
));

# Count data + meta data:
countTable <- read.csv(
  "results/tables/vsmcCounts.csv",
  header = TRUE
);

countTableDESeq <- sapply(countTable[, c(8:40)], as.integer);
row.names(countTableDESeq) <- countTable[, 1];
metadata <- data.frame(
  row.names = colnames(countTableDESeq),
  condition = condition,
  type = type,
  patient = patient
);

# Making DESeqDataSet object - adding design (correting for "patient"):
dds <- DESeqDataSetFromMatrix(
  countData = countTableDESeq,
  colData = metadata,
  design = ~ patient + condition
);

# Setting FBS02 as ref (can't use "contrast" due to running apeglm shrinkage):
#dds$condition <- relevel(dds$condition);
# Running DESeq2 with a Wald test:
dds <- DESeq(dds, test = "Wald");
rld <- rlog(dds);
assay(rld) <- limma::removeBatchEffect(assay(rld), batch = rld$patient, design = model.matrix(~condition, colData(rld)));
pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData = TRUE, ntop = 500);
names = rownames(colData(rld));
percentVar <- round(100 * attr(pcaData, "percentVar"));

pcaData$condition <- factor(
    metadata$condition, 
    levels=c(
        "Quiesced",
        "IL-1A/PDGF-BB",
        "Mock",
        "miR-CTRL",
        "miR-323",
        "miR-449b",
        "miR-491",
        "miR-892b",
        "miR-1827",
        "miR-4774",
        "miR-5681b"
    )
);

pdf(file = "results/figures/vsmc-pca.pdf")
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = patient)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#000000", "#B81059","#F03084", "#EA70A6","#224F9F","#426BB6","#678DD1","#81A8EF","#A4C3FB","#C7DAFC","#E0E8F8"))+
  ylim(-25,25)+ xlim(-35,15)+
  coord_fixed()+
  theme_1()
dev.off()

sampleDists<-dist(t(assay(rld)));
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-paste(rld$condition)
colnames(sampleDistMatrix)<-NULL
colors <- colorRampPalette( rev(brewer.pal(8, "BuPu")) )(255)
pdf(file = "results/figures/vsmc-sample-distances.pdf", width = 7, height = 6)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    fontsize_row = 9,
    show_rownames = TRUE,
    show_colnames = TRUE
)
dev.off()


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