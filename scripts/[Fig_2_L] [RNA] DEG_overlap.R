setwd("H:/Project/Screening Aging Genes/")
library(VennDiagram)

DEG <- read.table("./data/RNA-seq/20180613/DEG/P4_KR_NR.DEG.xls", header = T)
DEG <- DEG[DEG$logFC <= -0.585 & DEG$FDR <= 0.05, ]
DEG <- na.omit(DEG)

science.DEG <- read.delim("./reference/A Werner syndrome stem cell model unveils heterochromatin alterations as a driver of human aging/RNA-seq/hMSC.WTandKO.DEG.xls")
science.DEG <- science.DEG[science.DEG$logFC <= -0.585 & science.DEG$FDR <= 0.05, ]
science.DEG <- na.omit(science.DEG)

overlap <- intersect(as.character(science.DEG$gene), as.character(DEG$gene))

# write.table(overlap, "KAT7_downDEG.overlap.Science_upDEG.15foldchange.xls", row.names = F, col.names = F, sep = "\t", quote = F)
# write.table(overlap, "KAT7_upDEG.overlap.Science_downDEG.15foldchange.xls", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(overlap, "KAT7_downDEG.overlap.Science_downDEG.15foldchange.xls", row.names = F, col.names = F, sep = "\t", quote = F)


V <- venn.diagram(list(KAT7=as.character(DEG$gene), WRN=as.character(science.DEG$gene)),
                  filename=NULL, lwd=1, lty=2, 
                  fill=c("darkorange","dodgerblue"))

pdf("KAT7_downDEG.overlap.Science_downDEG.15foldchange.pdf", width = 5, height = 5)
grid.draw(V)
dev.off()
