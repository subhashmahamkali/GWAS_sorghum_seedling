#ml bedtools
#bedtools intersect -a -b  -wa -wb > out.bed

#library(Ropt)
library(data.table)
a=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/all_traits.txt", header=T,data.table=F)
a[,1]=gsub("_filtered.csv","",a[,1])
colnames(a)[1]="Trait"
ch=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt", header=T,data.table=F)
d=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/All_GWAS_sorg.txt", header=T,data.table=F)[,-c(5:6)]
d=merge(d,a,by="Trait",all.x=T)
d=d[,c(2,3,4,1,5:6)]

#k = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/GWAS/combined_GWAS.bed")
#t <- 6
#sum(-log10(k[, 4]) > t)
m = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/166_genes.txt")

col=c("slateblue","darkgreen","violetred","gold2","skyblue")

#col=c("slateblue4","#9D0191","#00C49A","darkgreen","navy")

col=adjustcolor(col,alpha.f = 0.5)
dd=NULL
for (k in 1:10){sub=subset(d,d[,1]==k)
  sub[,2]=sub[,2]+ch[k,4]
  dd=rbind(dd,sub)}
dd$POS=dd$POS/1e6 #convert from bp to mb

dd$col1=NA
dd[dd$Category=="Developmental",7]=col[2]
dd[dd$Category=="Seed",7]=col[4]
dd[dd$Category=="Architecture",7]=col[1]
dd[dd$Category=="Panicle",7]=col[5]
###overlap with B2
d2=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/positive_selection/NR_LN_HN_ove.txt", header=F,data.table=F)
d3=NULL
for (k in 1:10){sub=subset(d2,d2[,1]==k)
  sub[,5]=sub[,5]+ch[k,4]
  sub[,13]=sub[,13]+ch[k,4]
  sub[,21]=sub[,21]+ch[k,4]
  d3=rbind(d3,sub)
}
d3[,5]=d3[,5]/1e6 #column 5 and 13
d3[,13]=d3[,13]/1e6
d3[,21]=d3[,21]/1e6

###thr 282  5.928578
####hn
png("GWAS_HN_LN_NR_ov.png", height = 7, width = 14, res = 600, units = "in")
#pdf("sorg_GWAS_HN_LN_ov.pdf", height = 6, width = 10) 
#tiff("sorg_GWAS_HN_LN_ov.tiff",height = 120,width =130,res=600,units = "mm")
#par(mar=c(3.5,4,1,2), mfrow=c(2,1))
#par(mar=c(3,3,2,1), mfrow=c(3,1), oma=c(1,1,1,1))
par(mar = c(2, 3.5, 2, 2), mfrow = c(3,1), oma = c(2, 1, 1, 1))
##
d1=dd[dd[,5]=="HN",] #subset the HN from column 5
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,717),ylim=c(2,12),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
#this is used to plot manhattan plot
#d1[,2] used POS column for x-axis
#-log10(d1[,3]) used p values and converts to -log(10) and used as y-axis
#col=d1$col1 use the coresponding color to each dot
#pch=16 use solid circle markers for each SNP dot
#cex=0.4 represents the size

#bty="l",xlim=c(0,717),ylim=c(2,15),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)"
#This mentions the style of the plot #removing top and right borders for a cleaner plot
#x -axis limitations, y-axis limitations
#hide default axis and axis would can be added manually
#lable sixe to 60%
#x and y - axis labeles

#head(d1)
axis(2,las=2,tck=-.03,cex.axis=1.5,font.axis=2)
#making the Y-axis 
#adding the tick marks
#tck=-.03 - controlling the the tick mark size "+" = inside the plot and "-" = outside the plot
#cex.axis=0.8 controlls the tick size to 80% of the original

#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
#abline(h=5,col="red",lty=2)
segments(x0=0, x1=720, y0=5, y1=5, col="red", lty=2, lwd=2)
#drawing the horizantal line at -log(10)=5 and red color dashed line

abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.3))
#adding the vertical lines grey color as chromosome boundries
#mtext("Genomic Associations for Four Trait Categories Across Nitrogen Conditions (HN, LN, NR)", side=3, line=2.5, font=2, cex=1)#, outer = T)
#box()
legend("topleft",c("Architecture","Panicle","Developmental","Seed"),pch=16,col=col[c(1,5,2,4)],bty = "n",ncol=4,cex=1.5)
legend("topleft",c("Common Significant Loci"),pch=16,col=col[c(3)],bty = "n",ncol=1,cex=1.5,inset = c(0, 0.15))
#legend(x=360, y=15, "HN", text.col="red", bty="n", cex=1.2)
#adding the legend to the plot
#pch=16 defines the symbol as circle
#bty = "n" no box for the legend for cleaner appearance
#ncol=4,cex=1 arrange the legend with one row and 4 columns and default text size
mtext("High Nitrogen", side=3, line=0.8, font=2, cex=1,adj = 1)#, col = "firebrick")
points(d3[,5],-log10(d3[,4]),col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
#points are added with the pos from d3[,5] as x-axis
#-log10(d3[,4]) transform the p-value to -log10 as y-axis 
#adjusts the color wih transperancy of 0.7 and uses pch=8 star symbol and 0.8 size
#type="h" vertical lines drawn from baseline to -log10(p-value)
#lwd=1.3 increase the thickness for visibility


d1=dd[dd[,5]=="LN",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,717),ylim=c(2,12),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
mtext("-log10(p)", side=2, line=2.2, font=2, cex=1.2)

axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
#abline(h=5,col="red",lty=2)
segments(x0=0, x1=720, y0=5, y1=5, col="red", lty=2, lwd=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
mtext("Low Nitrogen", side=3, line=0.4, font=2, cex=1,adj = 1)#, col = "chartreuse4")
points(d3[,13],-log10(d3[,12]),col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)

d1=dd[dd[,5]=="NR",]


plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,720),ylim=c(2,12),axes=F,cex.lab=0.6,xlab="",ylab="", font.lab=2)
axis(2,las=2,tck=-.03,cex.axis=1.5, font.axis=2)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
#abline(h=5,col="red",lty=2)
segments(x0=0, x1=720, y0=5, y1=5, col="red", lty=2, lwd=2)

abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
mtext("Nitrogen response", side=3, line=0.5, font=2, cex=1,adj = 1)#, col = "yellow3")
points(d3[,21],-log10(d3[,20]),col=adjustcolor("violetred",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=1.5, font.axis=2)
mtext("Chromosome", side=1, line=2.5, font=2, cex=1.2)
#adds x-axis  and chromosome numbers at center position on each chromosome. with ticks also
#will display only in this panel
dev.off()

