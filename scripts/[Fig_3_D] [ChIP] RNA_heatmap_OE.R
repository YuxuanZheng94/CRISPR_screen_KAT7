setwd("H:/Project/Screening Aging Genes/")
library(gplots)

expression_raw <- read.table("./data/RNA-seq/20180715/merge.dexseq_clean.RPKM.xls", header = T, row.names = 1)
expression_log <- expression_raw[, c("LUC_rep1","LUC_rep2","KAT7_rep1","KAT7_rep2")]
# expression_raw <- read.table("./data/RNA-seq/20180613/merge.dexseq_clean.RPKM.xls", header = T, row.names = 1)
# expression_log <- expression_raw[, c("NR_P4_WT1","NR_P4_WT2","KR_P4_KO1","KR_P4_KO2")]
expression_log <- log2(expression_log + 1)

candidate.gene <- c("NNMT","ADAMTS14","SERPINE1","TPM1","THBS1","ANGPTL4","AMIGO2","SUGCT","LAMC2","TAGLN","MAMDC2",
                    "ADAM19","LMO7","CCDC80","CDKN2B","LOX","KCNG1","NUAK1","PLOD2","ANKRD1","COL4A1")
mydata <- expression_log[candidate.gene, ]

# plot
palette.breaks <- seq(-1, 1, 0.01)
# color.palette <- colorRampPalette(c("blue", "yellow"))(length(palette.breaks)-1)
color.palette <- colorRampPalette(c("blue","white","red"))(length(palette.breaks)-1)

pdf("./data/ChIP-seq/21gene/RNA/21gene_expression_OE.pdf", width = 3.5, height = 5)
heatmap.2(as.matrix(mydata), trace = "none", scale = "row", dendrogram = "none", 
          # Rowv = as.dendrogram(hr.candidates), Colv = "none",
          Rowv = "none", Colv = "none",
          # colsep = nrow(cluster.info[cluster.info$aging == "Y", ]), sepwidth = 2, sepcolor = "white",
          col = color.palette, breaks = palette.breaks, density.info = "none",
          cexCol = 1, cexRow = 0.6, labRow = NULL, key = T)
dev.off()
