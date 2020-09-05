setwd("H:/Project/Screening Aging Genes/")
library(ggplot2)

day28 <- read.table("./data/CRISPR-Cas9/D28.count/D28_vs_control.gene_summary_with_random_num.txt", header = T, row.names = 1)
day42 <- read.table("./data/CRISPR-Cas9/D42.count/D42_vs_control.gene_summary_with_random_num.txt", header = T, row.names = 1)
day56 <- read.table("./data/CRISPR-Cas9/D56.count/D56_vs_control.gene_summary_with_random_num.txt", header = T, row.names = 1)

plotfunction <- function(mydata)
{
  mydata <- mydata[order(mydata$pos.rank), ]

  mydata[1:20, "highlight"] <- 2
  mydata[is.na(mydata$highlight), "highlight"] <- 0
  mydata <- mydata[order(mydata$highlight), ]
  
  mydata[mydata$highlight == 1, "label"] <- rownames(mydata[mydata$highlight == 1, ])
  mydata[mydata$highlight == 2, "label"] <- rownames(mydata[mydata$highlight == 2, ])
  
  p <- ggplot(mydata, aes(x = gene.index, y = -log10(pos.p.value))) +
    # geom_point(aes(color = factor(highlight), size = factor(highlight))) +
    geom_point(aes(color = factor(highlight)), size = 1) +
    xlab("Gene index") + ylab("-log10 (P value)") +
    # scale_size_manual(values = c("0"=1, "1"=2, "2"=2), name = "") +
    scale_color_manual(values = c("0"="gray84","1"="darkslateblue","2"="red"), name = "") +
    geom_hline(yintercept = -log10(0.05), size = 1, color = "red") +
    # geom_text(aes(label=label), hjust = 0, vjust = 0, size = 1) +
    theme_classic()
  print(p)
}

pdf("CRISPR_gene_labeled_top20_gene_no_name.pdf", width = 5, height = 3)
plotfunction(mydata = day28)
plotfunction(mydata = day42)
plotfunction(mydata = day56)
dev.off()
