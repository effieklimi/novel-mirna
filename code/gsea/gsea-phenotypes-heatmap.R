library("ComplexHeatmap")
library("tidyverse")

options(warn = 1)
setwd("/Users/effieklimi/Documents/novel-mirna/")

draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)
}

# -------------------------------------- 
miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")

bp <- readRDS("results/rds/pathways/pathways-all-bioprocess.rds")
ke <- readRDS("results/rds/pathways/pathways-all-kegg.rds")
re <- readRDS("results/rds/pathways/pathways-all-reactome.rds")

# cell cycle
fgseaResCellCycle <- read.csv(
  "results/tables/fgsea-vsmc-cellcycle-sigNES-p01.csv",
  header = TRUE,
  row.names = 1
)

fgseaResCellCycle <- fgseaResCellCycle[order(rowSums(is.na(fgseaResCellCycle))), ]

fgseaResCellCycle[is.na(fgseaResCellCycle)] <- 0
fgseaRankedCellCycle <- as.matrix(fgseaResCellCycle[order(rowSums(fgseaResCellCycle), decreasing = FALSE),])
#matchesCellCycle <- which(rownames(fgseaResCellCycle) %in% rownames(fgseaRankedCellCycle[1:10,]))
fgseaMatrixCellCycle <- as.matrix(fgseaResCellCycle)
#rownames(fgseaMatrix) <- rownames(fgseaResCellCycle)
fgseaMatrixCellCycle <- fgseaMatrixCellCycle %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))


assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
anot <- data.frame(row.names = rownames(fgseaResCellCycle), Database = fgseaMatrixCellCycle[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "Gene Ontology BP" = "grey90"))
ha <- rowAnnotation(foo = anno_mark(at = matchesCellCycle, labels = rownames(fgseaMatrixCellCycle)[matches]))
pdf(file = "results/figures/fgsea-vsmc-heatmap-cellcycle.pdf", width = 10, height = 10)
ComplexHeatmap::pheatmap(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)],
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "grey",
    annotation_row = anot,
    annotation_colors = ann_colors,
    row_names_max_width = unit(10, "cm"),
    #main = "Gene Set Enrichment Analysis",
    colorRampPalette(c("#005dc7", "white", "#d60047"))(300),
    right_annotation = ha,
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)] == 0, "\\", ""), nrow(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()

# Motility
fgseaResMotility <- read.csv(
  "results/tables/fgsea-vsmc-motility-sigNES-p01.csv",
  header = TRUE,
  row.names = 1
)
fgseaResMotility[is.na(fgseaResMotility)] <- 0
fgseaRankedMotility <- as.matrix(fgseaResMotility[order(rowSums(fgseaResMotility), decreasing = FALSE),])
#matchesMotility <- which(rownames(fgseaResMotility) %in% rownames(fgseaRankedCellCycle[1:10,]))
fgseaMatrixMotility <- as.matrix(fgseaResMotility)
#rownames(fgseaMatrix) <- rownames(fgseaResMotility)
fgseaMatrixMotility <- fgseaMatrixMotility %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))


assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
anot <- data.frame(row.names = rownames(fgseaResMotility), Database = fgseaMatrixCellCycle[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "Gene Ontology BP" = "grey90"))
#ha <- rowAnnotation(foo = anno_mark(at = matchesCellCycle, labels = rownames(fgseaMatrixCellCycle)[matches]))
pdf(file = "results/figures/fgsea-vsmc-heatmap-motility.pdf")
ComplexHeatmap::pheatmap(fgseaMatrixMotility[,c(2,3,5,7,1,4,6)],
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "grey",
    annotation_row = anot,
    annotation_colors = ann_colors,
    #main = "Gene Set Enrichment Analysis",
    colorRampPalette(c("#005dc7", "white", "#d60047"))(300),
    #right_annotation = ha,
    border_gp = gpar(col = "black", lwd = 2),
    #display_numbers = matrix(ifelse(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)] == 0, "\\", ""), nrow(fgseaMatrixCellCycle[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()

# Deleterious phenotypes
fgseaResDelet <- read.csv(
  "results/tables/fgsea-vsmc-deleterious-sigNES-p01.csv",
  header = TRUE,
  row.names = 1
)
fgseaResDelet[is.na(fgseaResDelet)] <- 0
fgseaRankedDelet <- as.matrix(fgseaResDelet[order(rowSums(fgseaResDelet), decreasing = FALSE),])
matchesDelet <- which(rownames(fgseaResDelet) %in% rownames(fgseaRankedDelet[1:2,]))
fgseaMatrixDelet <- as.matrix(fgseaResDelet)
#rownames(fgseaMatrixDelet) <- rownames(fgseaResDelet)
fgseaMatrixDelet <- fgseaMatrixDelet %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  mutate(Database = case_when(
    as.character(pathway) %in% tolower(names(ke)) ~ "KEGG",
    as.character(pathway) %in% tolower(names(re)) ~ "Reactome",
    as.character(pathway) %in% tolower(names(bp)) ~ "Gene Ontology BP"
  )) %>%
  column_to_rownames(var = "pathway") %>%
  `colnames<-`(c("miR-1827", "miR-323a-3p", "miR-449b-5p", "miR-4774-3p", "miR-491-3p", "miR-5681b", "miR-892b", "Database"))



assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
anot <- data.frame(row.names = rownames(fgseaResultsSig), Database = fgseaMatrix[,8])
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "grey70", "Gene Ontology BP" = "grey90"))
ha <- rowAnnotation(foo = anno_mark(at = matches, labels = rownames(fgseaMatrix)[matches]))

pdf(file = "results/figures/fgsea-vsmc-heatmap-p01-multimir.pdf", width = 11, height = 10)
ComplexHeatmap::pheatmap(fgseaMatrix[,c(2,3,5,7,1,4,6)],
    border_color = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 7,
    fontsize_col = 7,
    show_rownames = FALSE,
    show_colnames = TRUE,
    na_col = "grey",
    annotation_row = anot,
    annotation_colors = ann_colors,
    main = "Gene Set Enrichment Analysis",
    colorRampPalette(c("#005dc7", "white", "#d60047"))(300),
    right_annotation = ha,
    border_gp = gpar(col = "black", lwd = 2),
    display_numbers = matrix(ifelse(fgseaMatrix[,c(2,3,5,7,1,4,6)] == 0, "\\", ""), nrow(fgseaMatrix[,c(2,3,5,7,1,4,6)]))
  )
  dev.off()