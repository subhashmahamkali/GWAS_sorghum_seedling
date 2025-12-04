library(data.table)
library(readxl)

setwd("~/Documents/gwas_sap")

#ch=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt", header=T,data.table=F)
ch=fread("data/2.2_positive_selec/sorg.chr_length_V5.txt", header=T,data.table=F)
#wild_vs_landrace
#d1=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg_wild_vs_landrace.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
d1=fread("data/2.2_positive_selec/sorg_wild_vs_landrace.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
d1$POS = (d1[,2] + d1[,3]) / 2
d1$POS = ifelse(d1$POS %% 1 == 0.5, floor(d1$POS), d1$POS)
d1 = d1[, -c(2,3)]
d1=d1[,c(1,3,2)]
d1=d1[order(d1[,1],d1[,2]),]
dp1=NULL
for (k in 1:10){sub=subset(d1,d1[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
dp1=rbind(dp1,sub)}
dp1$POS=dp1$POS/1e6
col1=ifelse(dp1[,1]%%2==1,"#00000066","#BEBEBE99")
thr1=quantile(dp1[,3],0.99)
thr1
#99% 
#0.4303338

#landrace_vs_improved
#d2=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/landrace_vs_improved.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
d2=fread("data/2.2_positive_selec/landrace_vs_improved.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
d2$POS = (d2[,2] + d2[,3]) / 2
d2$POS = ifelse(d2$POS %% 1 == 0.5, floor(d2$POS), d2$POS)
d2 = d2[, -c(2,3)]
d2=d2[,c(1,3,2)]
d2=d2[order(d2[,1],d2[,2]),]
dp2=NULL
for (k in 1:10){sub=subset(d2,d2[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
dp2=rbind(dp2,sub)}
dp2$POS=dp2$POS/1e6
col2=ifelse(dp2[,1]%%2==1,"#00000066","#BEBEBE99")
thr2=quantile(dp2[,3],0.99)
thr2
#99% 
#0.2395224

#ann <- read_excel(
  #"/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/W_I_L_positive_selection_new.xlsx",
  #sheet = "Sheet4"
#) |> as.data.frame()

ann <- read_excel(
  "data/2.2_positive_selec/W_I_L_positive_selection_new.xlsx",
  sheet = "Sheet4"
) |> as.data.frame()

if (!"chr"   %in% names(ann)) names(ann)[names(ann)=="V1"] <- "chr"
if (!"start" %in% names(ann)) names(ann)[names(ann)=="V2"] <- "start"
if (!"end"   %in% names(ann)) names(ann)[names(ann)=="V3"] <- "end"
ann$chr    <- as.numeric(ann$chr)
ann$start  <- as.numeric(ann$start)
ann$end    <- as.numeric(ann$end)
ann$mid_bp <- floor((ann$start + ann$end) / 2)
chr_col <- names(ch)[1]   
off_col <- names(ch)[4]   
offset_bp <- ch[[off_col]][ match(ann$chr, ch[[chr_col]]) ]
if (anyNA(offset_bp)) stop("Some ann$chr values not found in `ch`.")
ann$pos_mb <- (ann$mid_bp + offset_bp) / 1e6

#for (k in 1:10) {
  #sub <- subset(gene, gene$chr == k)
  #sub$pos <- sub$pos + ch[k, 4]
  #gene[gene$chr == k, "pos"] <- sub$pos
#}
#gene$pos <- gene$pos / 1e6
#png("1.pos.png", height = 7, width = 14, res = 600, units = "in")

tiff("graphs/01_publication/2.bal_pos_taj_pi/positive_selection.tiff", height = 7, width = 14, res = 600, units = "in", compression = "lzw", type = "cairo", bg = "white")
par(mar = c(2, 3.5, 2, 2), mfrow = c(2,1), oma = c(2, 2, 1, 1)) #outer margin area -BLTR
plot(dp1[,2],dp1[,3],col=col1,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("Fst", side=2, line=3.5, font=2, cex=1.2)
mtext("Wildtype vs. Landrace", side=3, line=0.4, font=2, cex=1,adj=1)
axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=1.2) 
#abline(h=thr,col="red",lty=2)
segments(x0=0, x1=720, y0=thr1, y1=thr1, col="red", lty=2, lwd=2)

tb1 <- subset(ann, annotation == "TB1")[1, ]
abline(v = tb1$pos_mb, col = "black", lty = 3)
y_tb1 <- max(dp1[,3], na.rm = TRUE) * 0.97   # place dot near panel top
points(tb1$pos_mb, y_tb1, pch = 21, bg = "dodgerblue2", col = "black", cex = 1)
#text(tb1$pos_mb, max(dp1[,3], na.rm = TRUE) * 0.95, labels = "TB1", font = 3, cex = 1.2)

plot(dp2[,2],dp2[,3],col=col2,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("Fst", side=2, line=3.5, font=2, cex=1.2)
mtext("Landrace vs. Improved", side=3, line=0.4, font=2, cex=1,adj=1)
axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=1.2) 
#abline(h=thr,col="red",lty=2)
segments(x0=0, x1=720, y0=thr2, y1=thr2, col="red", lty=2, lwd=2)
sh1 <- subset(ann, annotation == "sh1")[1, ]
abline(v = sh1$pos_mb, col = "black", lty = 3)
y_sh1 <- max(dp2[,3], na.rm = TRUE) * 0.97   # place dot near panel top
points(sh1$pos_mb, y_sh1, pch = 21, bg = "dodgerblue2", col = "black", cex = 1)
#text(sh1$pos_mb, max(dp2[,3], na.rm = TRUE) * 0.95, labels = "sh1", font = 3, cex = 1.2)
mtext("Chromosome", side=1, line=2, font=2, cex=1.2)
dev.off()
