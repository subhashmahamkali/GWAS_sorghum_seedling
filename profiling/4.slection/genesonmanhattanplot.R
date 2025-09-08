library(data.table)
library(readxl)




ch=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt", header=T,data.table=F)

gene = read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/pos_bal_loci_gene_annot.xlsx", sheet = 5)
gene$phys_pos <- as.integer(gene$phys_pos)
gene <- gene[gene$comparison != "sorg_wild_vs_improved", ]

for (k in 1:10) {
  sub <- subset(gene, gene$chr == k)
  sub$phys_pos <- sub$phys_pos + ch[k, 4]
  gene[gene$chr == k, "phys_pos"] <- sub$phys_pos
}
gene$phys_pos <- gene$phys_pos / 1e6

d=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)
d=d[,c(1,2,4)] #keeps only column 1,2,4 chr,pos,score
thr=quantile(d[,3],0.99)#threshold of 1%
thr
#1111.447
thrd=d[order(d[,1],d[,2]),]#order by chr and pos
dd=NULL
for (k in 1:10){sub=subset(d,d[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
dd=rbind(dd,sub)}
#adjust the genomic positiosn for all chromosomes continuesly 
dd$physPos=dd$physPos/1e6
#converts basepairs to megabase
col=ifelse(dd[,1]%%2==1,"#00000066","#BEBEBE99")
#defines color for odd and even chromosomes odd(black and even(grey)

png("/Users/subhashmahamkali/Documents/gwas_sap/1.b2_Fst_pos.png", height = 7, width = 14, res = 600, units = "in")
#pdf("b2_pos.pdf", height = 6, width = 10)  # Open PDF device
#library(Cairo)
#CairoSVG("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/b2_pos.svg", width = 10, height = 6)

#tiff("sorg_GWAS_HN_LN_NR2.tiff",height = 150,width =190,res=600,units = "mm")
par(mar = c(2, 3.5, 2, 2), mfrow = c(3,1), oma = c(2, 2, 1, 1))
#B L T R
#outer margin area -Bottem, Left, Top, Right
plot(dd[,2],dd[,3],col=col,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("B2", side=2, line=4, font=2, cex=1)
mtext("Balancing selection: wild population", side=3, line=0.4, font=2, cex=1,adj = 1)
axis(2,las=2,tck=-.03,cex.axis=1.5,font.axis=2)
segments(x0=0, x1=720, y0=thr, y1=thr, col="red", lty=2, lwd=2)
legend("topleft",c("Common Significant Loci"),pch=16,col="violetred",bty = "n",ncol=1,cex=1.5)

wl=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/bal_wild_lan.txt", header=F,data.table=F)
wl$V10 = wl[,2] + 40000 #creates new column of +40kb of start position 
wl$V11 = (wl[,6] + wl[,7])/2#average of start and end to plot the line shared between b2 and Fst
wl = wl[,c(1:3,10,4:7,11,8:9)]
#loops through the chr 1 - 10 and subsets rows per chromosome
wla=NULL
for (k in 1:10){sub=subset(wl,wl[,1]==k)
sub[,4]=sub[,4]+ch[k,4] 
sub[,9]=sub[,9]+ch[k,4] 
wla=rbind(wla,sub)}
wla$V10 =wla$V10/1e6
wla$V11 =wla$V11/1e6

li=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/bal_lan_imp.txt", header=F,data.table=F)
li$V10 = li[,2] + 40000
li$V11 = (li[,6] + li[,7])/2
li = li[,c(1:3,10,4:7,11,8:9)]
lia=NULL
for (k in 1:10){sub=subset(li,li[,1]==k)
sub[,4]=sub[,4]+ch[k,4] #GWAS
sub[,9]=sub[,9]+ch[k,4] #B2
lia=rbind(lia,sub)}
lia$V10 =lia$V10/1e6
lia$V11 =lia$V11/1e6
points(wla$V10,wla$V4,col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
points(lia$V10,lia$V4,col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
red_points1 <- gene[gene$comparison == "bal", ]  
points(red_points1$phys_pos, red_points1$score, col="red", pch=16, cex=1)


#wild_vs_landrace
d1=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg_wild_vs_landrace.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
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
col=ifelse(dp1[,1]%%2==1,"#00000066","#BEBEBE99")
thr=quantile(dp1[,3],0.99)
thr=quantile(dp1[,3],0.995)

plot(dp1[,2],dp1[,3],col=col,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("Fst", side=2, line=3.5, font=2, cex=1)
mtext("Positive selection: Wildtype vs. Landrace", side=3, line=0.4, font=2, cex=1,adj=1)
axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
segments(x0=0, x1=720, y0=thr, y1=thr, col="red", lty=2, lwd=2)
points(wla$V11,wla$V8,col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
red_points2 <- gene[gene$comparison == "sorg_wild_vs_landrace", ]  
points(red_points2$phys_pos, red_points2$score, col="red", pch=16, cex=1)


#landrace_vs_improved
d2=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/landrace_vs_improved.windowed.weir.fst", header=T,data.table=F)[,-c(4,6)]
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
col=ifelse(dp2[,1]%%2==1,"#00000066","#BEBEBE99")
thr=quantile(dp2[,3],0.99)
thr=quantile(dp2[,3],0.995)

plot(dp2[,2],dp2[,3],col=col,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("Fst", side=2, line=3.5, font=2, cex=1)
mtext("Positive Selection: Landrace vs. Improved", side=3, line=0.4, font=2, cex=1,adj=1)
axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
segments(x0=0, x1=720, y0=thr, y1=thr, col="red", lty=2, lwd=2)
points(lia$V11,lia$V8,col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=1.5,font.axis=2)
red_points3 <- gene[gene$comparison == "landrace_vs_improved", ]  
points(red_points3$phys_pos, red_points3$score, col="red", pch=16, cex=1)
mtext("Chromosome", side=1, line=2, font=2, cex=1)
dev.off()





