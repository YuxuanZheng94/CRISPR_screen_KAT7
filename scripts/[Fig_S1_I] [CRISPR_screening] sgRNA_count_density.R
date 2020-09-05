setwd("H:/Project/Screening Aging Genes/")
library(ggplot2)
library(pheatmap)
library(gplots)

sequencing.count <- c(18000109,11254450,9847561,17591196)

day28 <- read.table("./data/CRISPR-Cas9/D28.count/all.count.txt", header = T, row.names = 1)
day42 <- read.table("./data/CRISPR-Cas9/D42.count/all.count.txt", header = T, row.names = 1)
day56 <- read.table("./data/CRISPR-Cas9/D56.count/all.count.txt", header = T, row.names = 1)

day0 <- day28[, c("Gene","control.1.clean.cut","control.1.clean.cut")]

# normalization
day0[, 3] <- day0[, 3] / sequencing.count[1]
day28[, 3] <- day28[, 3] / sequencing.count[2]
day42[, 3] <- day42[, 3] / sequencing.count[3]
day56[, 3] <- day56[, 3] / sequencing.count[4]

plotfunction <- function(day.data)
{
  p <- ggplot(day.data, aes(x = log2(day.data[, 3]))) + 
    geom_histogram(binwidth = 0.1, fill = "lightslateblue", alpha = 0.7) +
    # geom_point(stat = "bin", binwidth = 0.1 ) +
    # geom_vline(xintercept = 10, linetype = "dashed", size = 1) +
    xlab("log2 (sgRNA Read Count)") + ylab("Count") +
    theme_classic()
  
  print(p)
}

pdf("sgRNA.count.distribution.pdf", width = 4.5, height = 4)
plotfunction(day.data = day0)
plotfunction(day.data = day28)
plotfunction(day.data = day42)
plotfunction(day.data = day56)
dev.off()



day0$day <- "day0"
day28$day <- "day28"
day42$day <- "day42"
day56$day <- "day56"

filter.day0 <- day0[, c(3,4)]; colnames(filter.day0) <- c("count","day")
filter.day28 <- day28[, c(3,4)]; colnames(filter.day28) <- c("count","day")
filter.day42 <- day42[, c(3,4)]; colnames(filter.day42) <- c("count","day")
filter.day56 <- day56[, c(3,4)]; colnames(filter.day56) <- c("count","day")

mydata <- rbind(filter.day0, filter.day28, filter.day42, filter.day56)
mydata <- mydata[mydata$count != 0, ]

a <- c("day10"="#00FFFF", "day28"="#C0FF3E", "day42"="#FF0000", "day56"="#1E90FF", "day0"="#FFB6C1")

pdf("./sgRNA_normalization_count_density_plot.pdf", width = 6, height = 5)
ggplot(mydata, aes(x = log2(count), color = day, fill = day)) + 
  scale_color_manual(values = a, name = "Stage") +
  scale_fill_manual(values = a, name = "Stage") +
  xlab("log2 (normalization count)") +
  geom_density(alpha = 0.2, size = 1) + theme_classic()

ggplot(mydata, aes(x = log2(count), color = day, fill = day)) + 
  scale_color_manual(values = a, name = "Stage") +
  scale_fill_manual(values = a, name = "Stage") +
  xlim(-10, 0) +
  xlab("log2 (normalization count)") +
  geom_density(alpha = 0.2, size = 1) + theme_classic()
dev.off()
