library(org.Mm.eg.db)
retrieved <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0030054", columns="ENSEMBL")