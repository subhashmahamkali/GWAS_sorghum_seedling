library(data.table)
d = fread("/Users/subhashmahamkali/Downloads/merged_chr1_maf.tsv")

plot(d$V3, d$V4,
     main = "Correlation between Zero-shot Score and MAF",
     xlab = "Zero-shot Score (deleteriousness)",
     ylab = "Minor Allele Frequency (MAF)",
     pch = 20, col = rgb(0.2, 0.4, 0.6, 0.4))
#lines(lowess(d$V3, d$V4), col = "red", lwd = 2)

d_sorted <- d[order(V3)]
top50 <- d_sorted[1:(.N/2)]
top10 <- d_sorted[1:floor(.N * 0.10)]
top1 <- d_sorted[1:floor(.N * 0.01)]
top0.1 <- d_sorted[1:floor(.N * 0.001)]
mean(top1$V4)
mean(top10$V4)
mean(top50$V4)
mean(top0.1$V4)
x_labels <- c("Top 50%", "Top 10%", "Top 1%", "Top 0.1%")
mean_maf <- c(mean(top50$V4), mean(top10$V4), mean(top1$V4), mean(top0.1$V4))
plot(mean_maf, type = "o", col = "steelblue", lwd = 2, pch = 16,
     xaxt = "n", ylim = c(0, 0.1),
     xlab = "Top deleterious variants)", ylab = "Mean MAF",
     main = "Mean MAF across increasing deleteriousness")
axis(1, at = 1:4, labels = x_labels)







library(ggplot2)
library(dplyr)


a = fread("/Users/subhashmahamkali/Downloads/snps_fixed_sorted.bed")
b = fread("/Users/subhashmahamkali/Downloads/all_features.bed")
unique(b$V4)


d = d[,-c (4)]
d$start = d$V2+1
d$end = d$V2-1
colnames(d) = c("chr","POS","zs","start","end")

# Read the file: assumes columns include chr, start, end, pos, score, ..., feature
#df <- read.table("/Users/subhashmahamkali/Downloads/annotated_snps_all.tsv", sep = "\t", header = FALSE)
#colnames(df) <- c("chr", "start", "end", "pos", "score", "chr_b", "start_b", "end_b", "feature", "dot", "strand")

# Boxplot
ggplot(df, aes(x = feature, y = score)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw() +
  labs(title = "Zero-Shot Deleteriousness by Feature Type",
       x = "Genomic Feature", y = "Zero-Shot Score")

