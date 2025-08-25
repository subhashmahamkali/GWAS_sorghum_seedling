#install.packages("R.utils")
library(data.table)
library(R.utils)

d=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/all_chr_merged_with_chr.B2.txt.gz", header=T,data.table=F)
d=d[,c(7,1,3)]
ch=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt", header=T,data.table=F)
d[,1]  <- as.numeric(d[,1])   # chromosome column
d[,2]  <- as.numeric(d[,2])   # position column
d[,3]  <- as.numeric(d[,3])
ch[,4] <- as.numeric(ch[,4])
dd=NULL
for (k in 1:10){sub=subset(d,d[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
dd=rbind(dd,sub)}
dd$physPos <- dd[,2] / 1e6 
midpoints <- sapply(1:10, function(k) {
  mean(range(dd[dd[,1] == k, 2], na.rm=TRUE)) / 1e6
})
dd_filt <- dd[dd[,3] >= 0, ]
thr <- quantile(dd_filt[,3], 0.99, na.rm = TRUE)
thr
#831.5203

la=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)
la=la[,c(1,2,4)] 
#thr <- quantile(la[,3], 0.99, na.rm = TRUE)
#thr
#1111.447
lan=NULL
for (k in 1:10){sub=subset(la,la[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
lan=rbind(lan,sub)}
lan$physPos <- lan[,2] / 1e6 
midpoins_l <- sapply(1:10, function(k) {
  mean(range(lan[lan[,1] == k, 2], na.rm=TRUE)) / 1e6
})
la_filt <- lan[lan[,3] >= 0, ]
thr_l <- quantile(la_filt[,3], 0.99, na.rm = TRUE)
thr_l
#1111.447

im=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/imp_chr_merged_with_chr.B2.txt.gz", header=T,data.table=F)
im=im[,c(7,1,3)]
im[,3] <- as.numeric(im[,3])
#thr <- quantile(im[,3], 0.99, na.rm = TRUE)
#thr
#1682.501
#ch[,4]  <- as.numeric(ch[,4])
im[,1]  <- as.numeric(im[,1])   # chromosome column
im[,2]  <- as.numeric(im[,2])   # position column
im[,3]  <- as.numeric(im[,3])
iim=NULL
for (k in 1:10){sub=subset(im,im[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
iim=rbind(iim,sub)}
#adjust the genomic positiosn for all chromosomes continuesly 
iim$physPos <- iim[,2] / 1e6 
#converts basepairs to megabase
midpoins_i <- sapply(1:10, function(k) {
  mean(range(iim[iim[,1] == k, 2], na.rm=TRUE)) / 1e6
})
ii_filt <- iim[iim[,3] >= 0, ]
thr_i <- quantile(ii_filt[,3], 0.99, na.rm = TRUE)
thr_i
#1682.746

c = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/W_L_I_overlap.bed")
w = c[,c(1:4)]
w$V5 <- (w$V2 + w$V3)/2 
w_plot <- NULL
for (k in 1:10) {
  sub = subset(w, w[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V5 = sub$V5 + ch[k,4]
  w_plot = rbind(w_plot, sub)
}
w_plot$V5 = w_plot$V5/1e6 

l = c[,c(5:8)]
l$V9 <- (l$V6 + l$V7)/2 
l_plot <- NULL
for (k in 1:10) {
  sub = subset(l, l[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V9 = sub$V9 + ch[k,4]
  l_plot = rbind(l_plot, sub)
}
l_plot$V9 = l_plot$V9/1e6 

i = c[,c(9:12)]
i$V13 <- (i$V10 + i$V11)/2 
i_plot <- NULL
for (k in 1:10) {
  sub = subset(i, i[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V13 = sub$V13 + ch[k,4]
  i_plot = rbind(i_plot, sub)
}
i_plot$V13 = i_plot$V13/1e6 


col_filt <- ifelse(dd_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")

#png("/Users/subhashmahamkali/Documents/gwas_sap/1.b2_pos.png", height = 7, width = 14, res = 600, units = "in")
png("main_plot.png",height = 7, width = 14, res = 600, units = "in")
legend("topleft",c("Common Significant Loci"),pch=16,col="violetred",bty = "n",ncol=1,cex=1.5)
par(mar = c(2, 3.5, 2, 2), mfrow = c(3,1), oma = c(2, 2, 1, 1))
plot(dd_filt$physPos, dd_filt[,3],
     col = col_filt, pch = 16, cex = 0.4,
     bty = "l", xlim = range(dd_filt$physPos, na.rm=TRUE),
     axes = FALSE, xlab = "", ylab = "", font.lab = 2)
# y-axis
axis(2, las=2, tck=-.03, cex.axis=1, font.axis=2)
# x-axis
chr_breaks <- seq(0, 717, length.out=11)
midpoints <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
axis(1, at=midpoints, labels=1:10, cex.axis=1.2, font.axis=2)
# threshold line
segments(x0=0, x1=720, y0=thr, y1=thr, col="red", lty=2, lwd=2)

#mtext("B2", side=2, line=3, font=2, cex=1.5)
#mtext("Chromosome", side=1, line=2.8, font=2, cex=1.5)
mtext("Balancing selection: wild population", side=3, line=0.5, font=2, cex=1, adj=1)
points(w_plot$V5, w_plot$V4, 
       col=adjustcolor("violetred", alpha.f=0.7),
       pch=8, cex=0.8, type="h", lwd=1.3)


col_la <- ifelse(la_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
plot(la_filt[,2], la_filt[,3],
     col = col_la,      # use the vector you just created
     pch = 16, cex = 0.4,
     bty = "l", xlim = c(0,717),
     axes = FALSE, xlab = "", ylab = "", font.lab = 2)
# y-axis
axis(2, las=2, tck=-.03, cex.axis=1, font.axis=2)
# x-axis: chromosomes
chr_breaks <- seq(0, 717, length.out=11)
midpoins_l <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
axis(1, at=midpoins_l, labels=1:10, cex.axis=1.2, font.axis=2)
segments(x0=0, x1=720, y0=thr_l, y1=thr_l, col="red", lty=2, lwd=2)
mtext("B2", side=2, line=3, font=2, cex=1.5)
#mtext("Chromosome", side=1, line=2.8, font=2, cex=1.5)
mtext("Balancing selection:landraces population", side=3, line=0.5, font=2, cex=1, adj=1)
points(l_plot$V9, l_plot$V8, 
       col=adjustcolor("violetred", alpha.f=0.7),
       pch=8, cex=0.8, type="h", lwd=1.3)



col_ii <- ifelse(ii_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
plot(ii_filt$physPos, ii_filt[,3], col = col_ii, pch = 16, cex = 0.4, bty = "l", xlim = c(0,717),
     axes = FALSE, xlab = "", ylab = "", font.lab = 2)
axis(2, las=2, tck=-.03, cex.axis=1, font.axis=2)
chr_breaks <- seq(0, 717, length.out=11)
midpoins_i <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
axis(1, at=midpoins_i, labels=1:10, cex.axis=1.2, font.axis=2)
segments(x0=0, x1=720, y0=thr_i, y1=thr_i, col="red", lty=2, lwd=2)
#mtext("B2", side=2, line=3, font=2, cex=1.5)
mtext("Chromosome", side=1, line=2.8, font=2, cex=1.5)
mtext("Balancing selection:improved population", side=3, line=0.5, font=2, cex=1, adj=1)
points(i_plot$V13, i_plot$V12, 
       col=adjustcolor("violetred", alpha.f=0.7),
       pch=8, cex=0.8, type="h", lwd=1.3)
dev.off() 



