library(Ropt)
library(data.table)
d=fread("propinqaum_3sample.frq", header=F,data.table=F)
a1=unlist(strsplit(d[,5],":"))
allele1=a1[seq(1,length(a1),by=2)]
freq1=as.numeric(a1[seq(2,length(a1),by=2)])
a2=unlist(strsplit(d[,6],":"))
allele2=a2[seq(1,length(a2),by=2)]
freq2=as.numeric(a2[seq(2,length(a2),by=2)])
d1=data.frame(d[,1:2],allele1,freq1,allele2,freq2)
d1[,3]=as.character(d1[,3])
d1[,5]=as.character(d1[,5])
an=apply(d1, 1, function(x){ifelse(x[4]>=x[6],x[3],x[5])})
ancestral_allele=an
d1=data.frame(d1,ancestral_allele)
colnames(d1)=c("Chr","Pos","allele1","Freq1","allele2","Freq2","ancestral_allele")
write.table(d1,file="Sorghum_ancestral_allele_V3.1.txt",row.names = F,quote = F,sep="\t",col.names = T)
