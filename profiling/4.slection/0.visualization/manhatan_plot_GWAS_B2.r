library(Ropt)
library(data.table)
a=fread("all_traits.txt", header=T,data.table=F)
a[,1]=gsub("_filtered.csv","",a[,1])
colnames(a)[1]="Trait"
ch=fread("sorg.chr_length_V5.txt", header=T,data.table=F)
d=fread("../sorghum_V5_B2/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)
d=d[,c(1,2,4)]
thr=quantile(d[,3],0.995)
d=d[order(d[,1],d[,2]),]
dd=NULL
for (k in 1:10)
{
  sub=subset(d,d[,1]==k)
  sub[,2]=sub[,2]+ch[k,4]
  dd=rbind(dd,sub)
}
dd$physPos=dd$physPos/1e6
col=ifelse(dd[,1]%%2==1,"#00000066","#BEBEBE99")
tiff("sorg_GWAS_HN_LN_NR2.tiff",height = 150,width =190,res=600,units = "mm")
par(mar=c(2,4,1,2),mfrow=c(5,1))

plot(dd[,2],dd[,3],col=col,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="B2")
axis(2,las=2,tck=-.03,cex.axis=0.8)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=thr,col="red",lty=2)

###overlap with B2
d2=fread("GWAS_sig_ov_B2.txt", header=T,data.table=F)
d3=NULL
for (k in 1:10)
{
  sub=subset(d2,d2[,1]==k)
  sub[,5]=sub[,5]+ch[k,4]
  sub[,11]=sub[,11]+ch[k,4]
  d3=rbind(d3,sub)
}
d3$MinP_pos=d3$MinP_pos/1e6
d3$B2_pos=d3$B2_pos/1e6
points(d3$B2_pos,d3[,13],col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)


d=fread("All_GWAS_sorg.txt", header=T,data.table=F)[,-c(5:6)]
d=merge(d,a,by="Trait",all.x=T)
d=d[,c(2,3,4,1,5:6)]
col=c("slateblue","darkgreen","violetred1","gold2","skyblue")
col=adjustcolor(col,alpha.f = 0.6)

dd=NULL
for (k in 1:10)
{
  sub=subset(d,d[,1]==k)
  sub[,2]=sub[,2]+ch[k,4]
  dd=rbind(dd,sub)
}
dd$POS=dd$POS/1e6
dd$col1=NA
dd[dd$Category=="Developmental",7]=col[2]
dd[dd$Category=="Seed",7]=col[4]
dd[dd$Category=="Architecture",7]=col[1]
dd[dd$Category=="Panicle",7]=col[5]


######hn
d1=dd[dd[,5]=="HN",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,717),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)")
axis(2,las=2,tck=-.03,cex.axis=0.8)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=5,col="red",lty=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
#box()
legend("topleft",c("Architecture","Panicle","Developmental","Seed"),pch=16,col=col[c(1,5,2,4)],bty = "n",ncol=4,cex=1)

hn1=d3[d3[,14]=="HN",]
points(hn1[,5],-log10(hn1[,4]),col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)

d1=dd[dd[,5]=="LN",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,720),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)")
axis(2,las=2,tck=-.03,cex.axis=0.8)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=5,col="red",lty=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
#box()
ln1=d3[d3[,14]=="LN",]
points(ln1[,5],-log10(ln1[,4]),col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)


d1=dd[dd[,5]=="NR",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,720),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)")
axis(2,las=2,tck=-.03,cex.axis=0.8)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=5,col="red",lty=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
#box()
nr1=d3[d3[,14]=="NR",]
points(nr1[,5],-log10(nr1[,4]),col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)

dev.off()
