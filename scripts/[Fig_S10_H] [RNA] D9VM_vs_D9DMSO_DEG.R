# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

setwd("H:/Project/Screening Aging Genes/")
library(DESeq2)
library(ggplot2)

expression.count <- read.table("./data/RNA-seq/20200824/merge.dexseq_clean.gene.xls", header = TRUE, row.names = 1)
cts <- expression.count[, grepl("D9_WM_", colnames(expression.count)) | grepl("D9_DMSO_", colnames(expression.count))]
colnames(cts)

coldata <- data.frame(array(NA, dim = c(ncol(cts), 3)))
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("Day","Group","Rep")
tmp <- strsplit(rownames(coldata), "_")
coldata$Day <- sapply(tmp, "[", 1)
coldata$Group <- sapply(tmp, "[", 2) 
coldata$Rep <- sapply(tmp, "[", 3) 
coldata$Group <- factor(coldata$Group, levels = c("DMSO","WM"))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Group)
dds

dds <- DESeq(dds)
res <- results(dds, pAdjustMethod = "BH")
# res

output <- data.frame(logFC = res$log2FoldChange, Pvalue = res$pvalue, FDR = res$padj, baseMean = res$baseMean, 
                     lfcSE = res$lfcSE, stat = res$stat, gene = rownames(cts))
mydata <- na.omit(output)

mydata[mydata$logFC <= -0.585 & mydata$FDR <= 0.05, "cluster"] <- "D9DMSO"
mydata[mydata$logFC >= 0.585 & mydata$FDR <= 0.05, "cluster"] <- "D9WM"
mydata[is.na(mydata$cluster), "cluster"] <- "pother"
mydata <- mydata[order(mydata$cluster, decreasing = T), ]
table(mydata$cluster)

pdf("D9DMSO_vs_D9WM_DEG_volcan.pdf", width = 6, height = 5)
ggplot(mydata, aes(x=mydata$logFC, y=-log10(mydata$FDR))) +
  geom_point(alpha = 0.8, aes(color = cluster, size = cluster)) +
  xlab("log2 (Fold Change)") + ylab("-log10 (FDR)") +
  scale_color_manual(values = c("D9DMSO"="blue","D9WM"="red","pother"="#D9D9D9")) +
  scale_size_manual(values = c("D9DMSO"=1.5,"D9WM"=1.5,"pother"=1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5)+
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", size = 0.5) +
  xlim(-6.4, 6.4) +
  theme_classic()
dev.off()

write.table(output, "D9DMSO_vs_D9WM_DEG.xls", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(mydata, "D9DMSO_vs_D9WM_DEG_clean.xls", row.names = F, col.names = T, sep = "\t", quote = F)



DEG_list <- mydata[mydata$cluster != "pother", ]
DEG_list <- DEG_list[order(DEG_list$cluster), ]

expression_RPKM <- read.table("./data/RNA-seq/20200824/merge.dexseq_clean.RPKM.xls", header = T, row.names = 1)
coldata <- coldata[order(coldata$Group), ]
heatmap_data <- expression_RPKM[as.character(DEG_list$gene), rownames(coldata)]

library(gplots)
palette.breaks <- seq(-1, 1, 0.01)
color.palette <- colorRampPalette(c("blue","white","red"))(length(palette.breaks)-1)

pdf("D9DMSO_vs_D9WM_DEG_heatmap.pdf", width = 5, height = 5)
heatmap.2(as.matrix(heatmap_data), trace = "none", scale = "row", dendrogram = "none",
          Rowv = FALSE, Colv = FALSE,
          # Rowv = as.dendrogram(hr.candidates), Colv = FALSE,
          # colsep = nrow(cluster.info[cluster.info$aging == "Y", ]), sepwidth = 2, sepcolor = "white",
          col = color.palette, breaks = palette.breaks, density.info = "none",
          cexCol = 1, cexRow = 0.6, labRow = NULL, key = T)
dev.off()
