setwd("H:/Project/Screening Aging Genes/")
library(gplots)

# select WRN DEG 
gene1 <- read.table("./data/RNA-seq/20200825/DEG/liver_2m_vs_28m_DEG_clean.xls", header = T)

gene1_up <- gene1[gene1$cluster == "28m", ]
gene1_down <- gene1[gene1$cluster == "2m", ]

# select KO DEG
gene2 <- read.table("./data/RNA-seq/20200825/DEG/liver_NTC_vs_KAT7_DEG_clean.xls", header = T)

gene2_up <- gene2[gene2$cluster == "KAT7", ]
gene2_down <- gene2[gene2$cluster == "NTC", ]

# define overlap
overlap.1 <- intersect(as.character(gene1_up$gene), as.character(gene2_down$gene))
overlap.2 <- intersect(as.character(gene1_down$gene), as.character(gene2_up$gene))

overlap.gene <- c(overlap.1, overlap.2)

# read WRN RNA-seq data
data1 <- read.table("./data/RNA-seq/20200825/merge.dexseq_clean.RPKM.xls", header = T, row.names = 1)
data2 <- read.table("./data/RNA-seq/20200825/merge.dexseq_clean.RPKM.xls", header = T, row.names = 1)

data1 <- data1[overlap.gene, grepl("Liver_2m_", colnames(data1)) | grepl("Liver_28m_", colnames(data1))]
data2 <- data2[overlap.gene, grepl("Liver_NTC_", colnames(data2)) | grepl("Liver_KAT7_", colnames(data2))]

data1 <- data1[, c("Liver_2m_1","Liver_2m_2","Liver_2m_3","Liver_28m_1","Liver_28m_2","Liver_28m_3")]
data2 <- data2[, c("Liver_NTC_1","Liver_NTC_2","Liver_NTC_3","Liver_KAT7_1","Liver_KAT7_2","Liver_KAT7_3")]

# read data into DESeq2
mydata <- cbind(data1, data2)
mydata <- log2(mydata + 1)

# plot
palette.breaks <- seq(-1, 1, 0.01)
color.palette <- colorRampPalette(c("blue","white","red"))(length(palette.breaks)-1)

pdf("heatmap_commone_DEG_between_2mVS28m_and_sgNTCVSsgKAT7.pdf", width = 3.5, height = 5)
heatmap.2(as.matrix(mydata), trace = "none", scale = "row", dendrogram = "none",
          Rowv = FALSE, Colv = FALSE,
          # Rowv = as.dendrogram(hr.candidates), Colv = FALSE,
          # colsep = nrow(cluster.info[cluster.info$aging == "Y", ]), sepwidth = 2, sepcolor = "white",
          col = color.palette, breaks = palette.breaks, density.info = "none",
          cexCol = 0.5, cexRow = 0.6, labRow = NULL, key = T)
dev.off()
