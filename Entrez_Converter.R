library("AnnotationDbi")
library("org.Hs.eg.db")

# RSV
data_RSV <- read.csv("../Network/Cytoscape/RSV_genes.csv", header = FALSE)

data_RSV$entrez = mapIds(org.Hs.eg.db,
                     keys=data_RSV$V1, #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

write(data_RSV$entrez,"../Network/Cytoscape/RSV_entrez.txt")

# HRV
data_HRV <- read.csv("../Network/Cytoscape/HRV_genes.csv", header = FALSE)

data_HRV$entrez = mapIds(org.Hs.eg.db,
                         keys=data_HRV$V1, #Column containing Ensembl gene ids
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

write(data_HRV$entrez,"../Network/Cytoscape/HRV_entrez.txt")

