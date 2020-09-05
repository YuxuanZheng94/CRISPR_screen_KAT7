# source("https://bioconductor.org/biocLite.R")
# biocLite("cummeRbund")
setwd("H:/Project/Screening Aging Genes/")
library(cummeRbund)

day28 <- read.table("./data/CRISPR-Cas9/D28.count/all.count.txt", header = T, row.names = 1)
day42 <- read.table("./data/CRISPR-Cas9/D42.count/all.count.txt", header = T, row.names = 1)
day56 <- read.table("./data/CRISPR-Cas9/D56.count/all.count.txt", header = T, row.names = 1)

abundance_function <- function(day.data)
{
  day.data$Gene <- as.character(day.data$Gene)
  day.data[grepl("NonTargetingControlGuideForHuman", day.data$Gene), "Gene"] <- "CONTROL"
  
  day.data_coltrol <- data.frame(sgRNA = rownames(day.data), gene = day.data[, 1], count = day.data[, 2])
  day.data_coltrol <- aggregate(day.data_coltrol$count, list(day.data_coltrol$gene), "sum")
  colnames(day.data_coltrol) <- c("gene","sum_count")
  day.data_coltrol$norm_count <- (day.data_coltrol$sum_count / sum(day.data[, 2])) * 10^6
  
  day.data_treat <- data.frame(sgRNA = rownames(day.data), gene = day.data[, 1], count = day.data[, 3])
  day.data_treat <- aggregate(day.data_treat$count, list(day.data_treat$gene), "sum")
  colnames(day.data_treat) <- c("gene","sum_count")
  day.data_treat$norm_count <- (day.data_treat$sum_count / sum(day.data[, 3])) * 10^6
  
  norm_day_data <- merge(day.data_coltrol, day.data_treat, by = "gene")
  colnames(norm_day_data) <- c("gene","control_count","control_norm","treat_count","treat_norm")
  
  norm_day_data$LFC <- log2(norm_day_data$treat_norm + 1) - log2(norm_day_data$control_norm + 1)
  # norm_day_data$LFC <- norm_day_data$treat_norm / (norm_day_data$control_norm + 0.000001)
  rownames(norm_day_data) <- norm_day_data$gene
  # norm_day_data <- norm_day_data[!grepl("NonTargetingControlGuideForHuman", rownames(norm_day_data)), ]

  return(norm_day_data)
}

FC_day28 <- abundance_function(day.data = day28)
FC_day42 <- abundance_function(day.data = day42)
FC_day56 <- abundance_function(day.data = day56)

# ########### tmp remove KAT7
# FC_day28 <- FC_day28[setdiff(rownames(FC_day28), "KAT7"), ]

FC_day42 <- FC_day42[rownames(FC_day28), ]
FC_day56 <- FC_day56[rownames(FC_day28), ]

new_data <- data.frame(D0 = 0, D28 = FC_day28$LFC, D42 = FC_day42$LFC, D56 = FC_day56$LFC)
rownames(new_data) <- rownames(FC_day28)

result <- genesCluster(object = new_data, k = 10, flag = "pos")
csClusterPlot(result, pseudocount = 0)
tmp_2 <- data.frame(cluster = result$clustering)
tmp_2$gene <- rownames(tmp_2)

tmp_3 <- tmp_2[grepl("Control", tmp_2$gene) | grepl("CONTROL", tmp_2$gene), ]
table(tmp_3$cluster)

write.table(tmp_2, "sgRNA_clustering_result_pos_thres2_cluster10.xls", row.names = F, col.names = T, sep = "\t", quote = F)
pdf("sgRNA_clustering_result_pos_thres2_cluster10.pdf", width = 10, height = 8)
csClusterPlot(result, pseudocount = 0)
dev.off()

genesCluster <- function(object, k, flag)
{
  require(cluster)
  m <- as.data.frame(object)

  # # exclude genes whose variance = 0
  # tmp <- data.frame(gene = rownames(m), var = apply(m[, c("D28","D42","D56")], 1, var))
  # m <- m[as.character(tmp[tmp$var != 0, "gene"]), ]
  
  if(flag == "neg")
  {
    # save genes who all less than 0
    m <- m[rowSums(m[, c("D28","D42","D56")] < 0) == 3, ]
  } else
  {
    # exclude genes who all less than 0
    m <- m[rowSums(m[, c("D28","D42","D56")] > 0) >= 1, ]
    
    tmp <- m[rowSums(m > 2) >= 1, ]
    m <- m[rownames(tmp), ]
  }
    
  # cluster
  n <- as.dist((1 - cor(t(m)))/2)
  clusters <- pam(n, k)
  class(clusters) <- "list"
  clusters$fpkm <- m
  
  return(clusters)
}


# day28 <- read.table("./data/CRISPR-Cas9/D28.count/D28_vs_control.sgrna_summary.txt", header = T, row.names = 1)
# day28 <- day28[!grepl("NonTargetingControlGuideForHuman", day28$Gene), ]
# day42 <- read.table("./data/CRISPR-Cas9/D42.count/D42_vs_control.sgrna_summary.txt", header = T, row.names = 1); day42 <- day42[rownames(day28), ]
# day56 <- read.table("./data/CRISPR-Cas9/D56.count/D56_vs_control.sgrna_summary.txt", header = T, row.names = 1); day56 <- day56[rownames(day28), ]
# 
# new_data <- data.frame(D28 = day28$LFC, D42 = day42$LFC, D56 = day56$LFC)
# rownames(new_data) <- rownames(day28)
# 
# depth_value <- colSums(new_data)
# 
# norm_new_data <- new_data
# for(j in 1:ncol(new_data))
# {
#   norm_new_data[, j] <- new_data[, j] / depth_value[j]
# }
# 
# write.table(data.frame(sgRNA = rownames(norm_new_data), norm_new_data), "normalization_sgRNA_abundance.xls", row.names = F, col.names = T, sep = "\t", quote = F)
# 
# norm_new_data <- read.table("/gpfs1/tangfuchou_pkuhpc/zhengyuxuan/RenJie/20190219/normalization_sgRNA_abundance.xls", row.names = 1, header = T)
# 
# tmp <- data.frame(sgRNA = rownames(norm_new_data), sum = rowSums(norm_new_data))
# for(i in 1:nrow(tmp))
# {
#   tmp[i, "max"] <- max(as.numeric(norm_new_data[i, ]))
# }
# tmp <- tmp[tmp$sum >= 0.00005, ]
# 
# # result <- genesCluster(object = norm_new_data[as.character(tmp$sgRNA), ],  k = 5)
# result <- genesCluster(object = new_data,  k = 10)
# csClusterPlot(result, pseudocount=0)
# tmp_2 <- data.frame(cluster = result$clustering)
# 
# tmp_2 <- cbind(tmp_2, day28[rownames(tmp_2), "Gene"])

library(ggplot2)

cluster_info <- read.table("./sgRNA_clustering_result_neg_cluster10.xls", header = T)
tmp_cluster_info <- cluster_info[cluster_info$cluster == 7, ]

tmp_mydata <- new_data[as.character(tmp_cluster_info$gene), ]

for(j in 1:ncol(tmp_mydata))
{
  tmp <- data.frame(gene = rownames(tmp_mydata), value = tmp_mydata[, j], time = colnames(tmp_mydata)[j])
  if(j == 1)
    mydata <- tmp else
      mydata <- rbind(mydata, tmp)
}

# mydata[mydata$value < 0, "value"] <- 0
mydata$color <- "gene"
# mydata[mydata$gene == "SUV39H1", "color"] <- "zhighlight"
# mydata$gene <- as.character(mydata$gene)
# mydata[mydata$gene == "SUV39H1", "gene"] <- "zSUV39H1"
mydata[mydata$gene == "CONTROL", "color"] <- "zhighlight"
mydata$gene <- as.character(mydata$gene)
mydata[mydata$gene == "CONTROL", "gene"] <- "zCONTROL"

tmp_mean <- aggregate(mydata$value, list(mydata$time), "mean")
tmp_mean <- data.frame(gene = "zmean", value = tmp_mean$x, time = tmp_mean$Group.1, color = "zmean")
mydata <- rbind(mydata, tmp_mean)

mydata <- mydata[order(mydata$color), ]

pdf("highlight_cluster7_CONTROL.pdf", width = 4, height = 3)
ggplot(mydata, aes(x = time, y = value, group = gene)) +
  geom_line(aes(color = color, size = color)) +
  scale_color_manual(values = c("gene"="cyan","zhighlight"="#6495ED","zmean"="black")) +
  scale_size_manual(values = c("gene"=0.5,"zhighlight"=2,"zmean"=1)) +
  theme_classic()
dev.off()

pdf("highlight_cluster10_KAT7.pdf", width = 4, height = 3)
ggplot(mydata, aes(x = time, y = value, group = gene)) +
  geom_line(aes(color = color, size = color)) +
  scale_color_manual(values = c("gene"="cyan","zhighlight"="#6495ED","zmean"="black")) +
  scale_size_manual(values = c("gene"=0.5,"zhighlight"=2,"zmean"=1)) +
  theme_classic()
dev.off()
