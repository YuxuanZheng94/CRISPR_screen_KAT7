# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

setwd("H:/Project/Screening Aging Genes/reference/A Werner syndrome stem cell model unveils heterochromatin alterations as a driver of human aging/RNA-seq/")
library(DESeq2)
library(ggplot2)

expression.count <- read.table("./merge.read.count.xls", header = TRUE, row.names = 1)
cts <- expression.count
colnames(cts)

coldata <- data.frame(array(NA, dim = c(ncol(cts), 3)))
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("Group","Type","Rep")
tmp <- strsplit(rownames(coldata), "_")
coldata$Group <- sapply(tmp, "[", 1)
coldata$Type <- sapply(tmp, "[", 2)
coldata$Rep <- rep(c(1,2)) 
coldata$Group <- factor(coldata$Group, levels = c("WT","KO"))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Group)
# dds

dds <- DESeq(dds)
res <- results(dds, pAdjustMethod = "BH")
# res

output <- data.frame(logFC = res$log2FoldChange, Pvalue = res$pvalue, FDR = res$padj, baseMean = res$baseMean, 
                     lfcSE = res$lfcSE, stat = res$stat, gene = rownames(cts))
mydata <- na.omit(output)

mydata[mydata$logFC <= -0.585 & mydata$FDR <= 0.05, "cluster"] <- "WT"
mydata[mydata$logFC >= 0.585 & mydata$FDR <= 0.05, "cluster"] <- "KO"
mydata[is.na(mydata$cluster), "cluster"] <- "zother"
mydata <- mydata[order(mydata$cluster, decreasing = T), ]

pdf("hMSC.WTandKO.DEG.pdf", width = 6, height = 5)
ggplot(mydata, aes(x=mydata$logFC, y=-log10(mydata$FDR))) +
  geom_point(alpha = 0.8, aes(color = cluster, size = cluster)) + 
  xlab("log2 (Fold Change)") + ylab("-log10 (FDR)") + 
  scale_color_manual(values = c("WT"="dodgerblue","KO"="darkorange","zother"="#D9D9D9")) +
  scale_size_manual(values = c("WT"=1.5,"KO"=1.5,"zother"=1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5)+ 
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", size = 0.5) +
  # xlim(-1.4, 1.4) +
  theme_classic()
dev.off()

write.table(output, "hMSC.WTandKO.DEG.xls", row.names = F, col.names = T, sep = "\t", quote = F)
