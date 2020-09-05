setwd("H:/Project/Screening Aging Genes/")
library(ggplot2)

mydata <- read.table("./data/RNA-seq/20180715/DEG/LUC_KAT7.DEG.xls", header = T)
mydata <- na.omit(mydata)

mydata[mydata$logFC <= -0.585 & mydata$FDR <= 0.05, "cluster"] <- "NR"
mydata[mydata$logFC >= 0.585 & mydata$FDR <= 0.05, "cluster"] <- "KR"
mydata[is.na(mydata$cluster), "cluster"] <- "pother"
mydata <- mydata[order(mydata$cluster, decreasing = T), ]

pdf("KAT7_OE_DEG_new_color.pdf", width = 6, height = 5)
ggplot(mydata, aes(x=mydata$logFC, y=-log10(mydata$FDR))) +
  geom_point(alpha = 0.8, aes(color = cluster, size = cluster)) +
  xlab("log2 (Fold Change)") + ylab("-log10 (FDR)") +
  scale_color_manual(values = c("NR"="blue","KR"="red","pother"="#D9D9D9")) +
  scale_size_manual(values = c("NR"=1,"KR"=1,"pother"=1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5)+
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", size = 0.5) +
  xlim(-3.5, 3.5) +
  theme_classic()
dev.off()
