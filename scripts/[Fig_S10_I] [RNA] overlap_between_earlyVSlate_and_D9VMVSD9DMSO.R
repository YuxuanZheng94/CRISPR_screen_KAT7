setwd("H:/Project/Screening Aging Genes/")
library(VennDiagram)

gene1 <- read.table("./data/RNA-seq/20200824/DEG/D4DMSO_vs_D9DMSO_DEG_clean.xls", header = T)
gene1_up <- gene1[gene1$cluster == "D9", ]
gene1_down <- gene1[gene1$cluster == "D4", ]

gene2 <- read.table("./data/RNA-seq/20200824/DEG/D9DMSO_vs_D9WM_DEG_clean.xls", header = T)
gene2_up <- gene2[gene2$cluster == "D9WM", ]
gene2_down <- gene2[gene2$cluster == "D9DMSO", ]

tmp1 <- intersect(as.character(gene1_up$gene), as.character(gene2_down$gene))
write.table(tmp1, "LateUP_and_WMDOWN.xls", row.names = F, col.names = F, sep = "\t", quote = F)
tmp2 <- intersect(as.character(gene1_down$gene), as.character(gene2_up$gene))
write.table(tmp2, "LateDOWN_and_WMUP.xls", row.names = F, col.names = F, sep = "\t", quote = F)

V <- venn.diagram(list(gene1_up=as.character(gene1_up$gene), gene2_down=as.character(gene2_down$gene)),
                  filename=NULL, lwd=1, lty=2,
                  fill=c("darkorange","dodgerblue"))
V <- venn.diagram(list(gene1_down=as.character(gene1_down$gene), gene2_up=as.character(gene2_up$gene)),
                  filename=NULL, lwd=1, lty=2,
                  fill=c("darkorange","dodgerblue"))

# pdf("LateUP_and_WMDOWN.pdf", width = 5, height = 5)
pdf("LateDOWN_and_WMUP.pdf", width = 5, height = 5)
grid.draw(V)
dev.off()
