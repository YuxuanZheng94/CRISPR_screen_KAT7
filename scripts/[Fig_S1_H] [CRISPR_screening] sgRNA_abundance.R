setwd("H:/Project/Screening Aging Genes/")
library(ggplot2)
library(pheatmap)
library(gplots)

day28 <- read.table("./data/CRISPR-Cas9/D28.count/all.count.txt", header = T, row.names = 1)
day42 <- read.table("./data/CRISPR-Cas9/D42.count/all.count.txt", header = T, row.names = 1)
day56 <- read.table("./data/CRISPR-Cas9/D56.count/all.count.txt", header = T, row.names = 1)
# day56merge <- read.table("./data/CRISPR-Cas9/D56merge.count/all.count.txt", header = T, row.names = 1)

day0 <- day28[, c("Gene","control.1.clean.cut","control.1.clean.cut")]

# get top rank genes of terminal day
# gene.rank <- read.table("./data/CRISPR-Cas9/D56merge_vs_control.gene_summary.txt", header = T)
gene.rank <- read.table("./data/CRISPR-Cas9/D56.count/D56_vs_control.gene_summary.txt", header = T)
gene.rank <- gene.rank[order(gene.rank$pos.rank), ]

top.gene <- as.character(gene.rank[1:16, "id"])
top.gene <- top.gene[!grepl("Control", top.gene)]

# calculate sum normalization sgRNA read count 
sumcount <- function(day.data)
{
  # calculate Top 15 sgRNA
  for(i in 1:length(top.gene))
  {
    tmp <- day.data[day.data$Gene == top.gene[i], ]
    tmp <- sum(tmp[, 3]) / sum(day.data[, 3])
    tmp <- data.frame(abundance = tmp, gene = top.gene[i])
    
    if(i == 1)
      mydata <- tmp else
        mydata <- rbind(mydata, tmp)
  }
  
  # calculate control sgRNA
  tmp <- day.data[grepl("Control", day.data$Gene), ]
  tmp <- sum(tmp[, 3]) / sum(day.data[, 3])
  tmp <- data.frame(abundance = tmp, gene = "control")
  mydata <- rbind(mydata, tmp)
  
  return(mydata)
}

abundance.d0 <- sumcount(day.data = day0)
abundance.d28 <- sumcount(day.data = day28)
abundance.d42 <- sumcount(day.data = day42)
abundance.d56 <- sumcount(day.data = day56)
# abundance.d56merge <- sumcount(day.data = day56merge)

mydata <- rbind(abundance.d0, abundance.d28, abundance.d42, abundance.d56)
mydata$day <- c(rep("day0",16), rep("day28",16), rep("day42",16), rep("day56",16))
# mydata <- rbind(abundance.d0, abundance.d28, abundance.d42, abundance.d56merge)
# mydata$day <- c(rep("day0",16), rep("day28",16), rep("day42",16), rep("day56merge",16))
mydata$info <- paste(mydata$day, mydata$gene, sep = "_")

# paint
pdf("Top15.sgRNA.abundance.pdf", width = 10, height = 4)
ggplot(mydata, aes(x = factor(info, levels = mydata$info), y = abundance)) + 
  geom_bar(position=position_dodge(), stat="identity", width = 0.7, fill = "dodgerblue")  + 
  xlab("") + ylab("Abundance of sgRNA") +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45))
dev.off()

write.table(mydata, "Top15.sgRNA.abundance.xls", row.names = F, col.names = T, sep = "\t", quote = F)




heatmap.mydata <- data.frame(array(NA, dim = c(16,4)))
rownames(heatmap.mydata) <- c(top.gene, "control")
colnames(heatmap.mydata) <- c("day0","day28","day42","day56")

for(i in 1:nrow(heatmap.mydata))
{
  for(j in 1:ncol(heatmap.mydata))
  {
    heatmap.mydata[i, j] <- mydata[mydata$gene == rownames(heatmap.mydata)[i] & mydata$day == colnames(heatmap.mydata)[j], "abundance"]
  }
}

pdf("./Top15_sgRNA_abundance_heatmap.pdf", width = 10, height = 10)
pheatmap(heatmap.mydata, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, color =
           # colorRampPalette(c("#FFFDD2", "#FFDF84","#FF9D3C","#FF0000","#B5002D"))(50),
           colorRampPalette(c("#2A318C","#F8EE18","#EE6F6E","#EE3B3B","#CD3333","#8B2323","#8E181B"))(30),
           # colorRampPalette(c("blue4","white","#C22625"))(50),
         border_color = "grey", cellwidth = 20, cellheight = 20)

pheatmap(heatmap.mydata, scale = "row", cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, color =
           # colorRampPalette(c("#FFFDD2", "#FFDF84","#FF9D3C","#FF0000","#B5002D"))(length(seq(-2, 2, 0.01))-1),
           colorRampPalette(c("#2A318C","#F8EE18","#EE6F6E","#EE3B3B","#CD3333","#8B2323","#8E181B"))(length(seq(-2, 2, 0.01))-1),
         breaks = seq(-2, 2, 0.01),
         border_color = "grey", cellwidth = 20, cellheight = 20)
dev.off()

