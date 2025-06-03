library(data.table)
#d=fread("/Users/subhashmahamkali/Downloads/output_2.tsv", header=T,data.table=F)
d=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/merged_output.tsv", header=T,data.table=F)


d=d[,c(1,4,8)]
thr=quantile(d[,3],0.99)
d=d[order(d[,1],d[,2]),]
d$pos=d$pos/1e6
filtered_data <- d[d$pos <= 30, ]

plot(filtered_data$pos, filtered_data$zeroShotScore,
     pch = 20, col = "darkblue",
     xlab = "chr1(Mb)",
     ylab = "Zero-Shot Score value",
     main = "deleterious mutation")
threshold <- quantile(d$zeroShotScore, 0.01)  # bottom 1%
abline(h = threshold, col = "red", lty = 2)
deleterious <- d$zeroShotScore < threshold
points(d$pos[deleterious], d$zeroShotScore[deleterious],
       col = "red", pch = 20)







negative_data <- filtered_data[filtered_data$zeroShotScore < 0,]

top_50_negative <- quantile(negative_data$zeroShotScore, 0.50)
top_10_negative <- quantile(negative_data$zeroShotScore, 0.10)
top_1_negative <- quantile(negative_data$zeroShotScore, 0.01)
top_01_negative <- quantile(negative_data$zeroShotScore, 0.001)

hist(filtered_data$zeroShotScore, 
     breaks = 50, 
     main = "Histogram of zeroShotScore (Positive and Negative)",
     xlab = "zeroShotScore",
     col = "lightblue", 
     border = "white", 
     xlim = c(min(filtered_data$zeroShotScore), max(filtered_data$zeroShotScore)),
     ylab = "Frequency")

abline(v = top_50_negative, col = "red", lwd = 2, lty = 2)
abline(v = top_10_negative, col = "green", lwd = 2, lty = 2)
abline(v = top_1_negative, col = "blue", lwd = 2, lty = 2)
abline(v = top_01_negative, col = "purple", lwd = 2, lty = 2)

text(top_50_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "50%", pos = 4, col = "red")
text(top_10_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "10%", pos = 4, col = "green")
text(top_1_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "1%", pos = 4, col = "blue")
text(top_01_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "0.1%", pos = 4, col = "purple")




















# Open a PNG device to save the plot
png("zeroShotScore_histogram.png", width = 800, height = 600, res = 00)

# Create the histogram
o = hist(filtered_data$zeroShotScore, 
         breaks = 50, 
         main = "Histogram of zeroShotScore (Positive and Negative)",
         xlab = "zeroShotScore",
         col = "gold", 
         border = "white", 
         xlim = c(min(filtered_data$zeroShotScore), max(filtered_data$zeroShotScore)),
         ylab = "Frequency")

# Add vertical lines for quantiles
abline(v = top_50_negative, col = "red", lwd = 2, lty = 2)
abline(v = top_10_negative, col = "green", lwd = 2, lty = 2)
abline(v = top_1_negative, col = "blue", lwd = 2, lty = 2)
abline(v = top_01_negative, col = "purple", lwd = 2, lty = 2)

# Add labels for the quantiles
text(top_50_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "50%", pos = 4, col = "red")
text(top_10_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "10%", pos = 4, col = "green")
text(top_1_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "1%", pos = 4, col = "blue")
text(top_01_negative, max(hist(filtered_data$zeroShotScore, breaks = 50, plot = FALSE)$counts), 
     labels = "0.1%", pos = 4, col = "purple")

# Close the PNG device (this saves the file)
dev.off()
