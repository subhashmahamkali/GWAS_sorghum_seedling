library(Ropt)
library(data.table)
d=fread("Sorghum_ancestral_allele_V3.1.txt", header=T,data.table=F)
for(i in 1:10)
{
  d1=d[d[,1]==i,]
  out=qq("Chr{i}.Sorghum_ancestral_allele_V3.1.txt")
  write.table(d1,file=out,row.names = F,quote = F,sep="\t",col.names = T)
  
}