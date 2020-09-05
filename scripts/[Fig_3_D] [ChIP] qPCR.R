setwd("H:/Project/Screening Aging Genes/data/ChIP-seq/21gene/qPCR/")
library(gplots)

expression_raw <- read.table("qPCR.txt", header = T, row.names = 1)

candidate.gene <- c("NNMT","ADAMTS14","SERPINE1","TPM1","THBS1","ANGPTL4","AMIGO2","SUGCT","LAMC2","TAGLN","MAMDC2",
                    "ADAM19","LMO7","CCDC80","CDKN2B","LOX","KCNG1","NUAK1","PLOD2","ANKRD1","COL4A1")
mydata <- expression_raw[candidate.gene, ]

# plot
palette.breaks <- seq(0, 1, 0.01)
# color.palette <- colorRampPalette(c("blue", "yellow"))(length(palette.breaks)-1)
color.palette <- colorRampPalette(c("blue","white","red"))(length(palette.breaks)-1)

pdf("21gene_qPCR.pdf", width = 3.5, height = 5)
heatmap.2(as.matrix(mydata), trace = "none", scale = "none", dendrogram = "none", 
          # Rowv = as.dendrogram(hr.candidates), Colv = "none",
          Rowv = "none", Colv = "none",
          # colsep = nrow(cluster.info[cluster.info$aging == "Y", ]), sepwidth = 2, sepcolor = "white",
          col = color.palette, breaks = palette.breaks, density.info = "none",
          cexCol = 1, cexRow = 0.6, labRow = NULL, key = T)
dev.off()
