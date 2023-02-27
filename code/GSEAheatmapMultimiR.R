##############################
library(pheatmap)
##############################

draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)
}


fgseaResultsSig50Top2 <- read.csv(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResultsSig50Top2NES.csv", 
  header = TRUE, row.names = 1)

fgseaResultsSig50 <- read.csv(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResultsSig50NES.csv", 
  header = TRUE, row.names = 1)

fgseaResultsSig100 <- read.csv(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/multimiR/GSEA multimiR/fgseaResultsSig100NES.csv", 
  header = TRUE, row.names = 1)


# Rank results based on decreasing number of NAs per row and format for heatmap:
fgseaResultsSig50Top2NES <- as.matrix(fgseaResultsSig50Top2[order(rowSums(is.na(fgseaResultsSig50Top2))), c(2, 3, 5, 7, 1, 4, 6)])
fgseaResultsSig50NES <- as.matrix(fgseaResultsSig50[order(rowSums(is.na(fgseaResultsSig50))), c(2, 3, 5, 7, 1, 4, 6)])
fgseaResultsSig100NES <- as.matrix(fgseaResultsSig100[order(rowSums(is.na(fgseaResultsSig100))), c(2, 3, 5, 7, 1, 4, 6)])


assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
ann_colors <- list(Database = c("KEGG" = "#636362", "Reactome" = "#babab8", "Gene Ontology BP" = "white"))

anot50Top2 <- data.frame(row.names = rownames(fgseaResultsSig50Top2), Database = fgseaResultsSig50Top2$Database)
anot50 <- data.frame(row.names = rownames(fgseaResultsSig50), Database = fgseaResultsSig50$Database)
anot100 <- data.frame(row.names = rownames(fgseaResultsSig100), Database = fgseaResultsSig100$Database)


pheatmap(fgseaResultsSig50Top2NES,
         border_color = "black",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row = 8,
         fontsize_col = 5,
         show_rownames = TRUE,
         show_colnames = TRUE,
         na_col = "white",
         annotation_row = anot50Top2,
         annotation_colors = ann_colors,
         main = "Gene Set Enrichment Analysis",
         colorRampPalette(c("#005dc7", "white", "#d60047"))(300))

pheatmap(fgseaResultsSig50NES,
         border_color = "black",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row = 7,
         fontsize_col = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         na_col = "white",
         annotation_row = anot50,
         annotation_colors = ann_colors,
         main = "Gene Set Enrichment Analysis",
         colorRampPalette(c("#005dc7", "white", "#d60047"))(300))

pheatmap(fgseaResultsSig100NES,
         border_color = "black",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row = 7,
         fontsize_col = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         na_col = "white",
         annotation_row = anot100,
         annotation_colors = ann_colors,
         main = "Gene Set Enrichment Analysis",
         colorRampPalette(c("#005dc7", "white", "#d60047"))(300))

