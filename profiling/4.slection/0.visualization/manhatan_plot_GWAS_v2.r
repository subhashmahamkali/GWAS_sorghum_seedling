library(Ropt)
library(data.table)
a=fread("all_traits.txt", header=T,data.table=F)
a[,1]=gsub("_filtered.csv","",a[,1])
colnames(a)[1]="Trait"
ch=fread("sorg.chr_length_V5.txt", header=T,data.table=F)
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
###overlap with B2
d2=fread("GWAS_HN_LN_overlaps.txt", header=F,data.table=F)
d3=NULL
for (k in 1:10)
{
  sub=subset(d2,d2[,1]==k)
  sub[,5]=sub[,5]+ch[k,4]
  sub[,13]=sub[,13]+ch[k,4]
  d3=rbind(d3,sub)
}
d3[,5]=d3[,5]/1e6
d3[,13]=d3[,13]/1e6

###thr 282  5.928578
####hn
tiff("sorg_GWAS_HN_LN_ov.tiff",height = 120,width =130,res=600,units = "mm")
par(mar=c(2,4,1,2),mfrow=c(3,1))
##
d1=dd[dd[,5]=="HN",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,717),ylim=c(2,15),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)")
axis(2,las=2,tck=-.03,cex.axis=0.8)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=5,col="red",lty=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
#box()
legend("topleft",c("Architecture","Panicle","Developmental","Seed"),pch=16,col=col[c(1,5,2,4)],bty = "n",ncol=4,cex=1)

points(d3[,5],-log10(d3[,4]),col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)

d1=dd[dd[,5]=="LN",]
plot(d1[,2],-log10(d1[,3]),col=d1$col1,pch=16,cex=0.4,bty="l",xlim=c(0,720),ylim=c(2,15),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)")
axis(2,las=2,tck=-.03,cex.axis=0.8)
#axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
abline(h=5,col="red",lty=2)
abline(v=ch[,4]/1e6,col=adjustcolor("gray",alpha.f = 0.6))
#box()
points(d3[,13],-log10(d3[,12]),col=adjustcolor("violetred1",alpha.f = 0.7),pch=8,cex=0.8,type="h",lwd=1.3)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=0.8) 
dev.off()
