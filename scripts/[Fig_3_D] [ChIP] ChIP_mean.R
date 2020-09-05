setwd("H:/Project/Screening Aging Genes/")

chip_data <- read.delim("./data/ChIP-seq/21gene_ChIP/H3K14ac_21genes_matrix_noheader.txt", header = F)
basic_info <- chip_data[, 1:6]
chip_data <- chip_data[, 7:ncol(chip_data)]

sample_1 <- chip_data[, 1:600] # NC_P4_H3K14ac_rep1
sample_2 <- chip_data[, 601:1200] # NC_P4_H3K14ac_rep3
sample_3 <- chip_data[, 1201:1800] # KC_P4_H3K14ac_rep2
sample_4 <- chip_data[, 1801:2400] # KC_P4_H3K14ac_rep3

new_sample_1 <- (sample_1 + sample_2) / 2 # mean of NC_H3K14ac
new_sample_2 <- (sample_3 + sample_4) / 2 # mean of KC_H3K14ac

new_chip_data <- cbind(basic_info, new_sample_1, new_sample_2)

write.table(new_chip_data, "./data/ChIP-seq/21gene_ChIP/H3K14ac_21genes_mean_matrix", row.names = F, col.names = F, sep = "\t", quote = F)
