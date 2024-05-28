library(lme4)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(Matrix)
library(data.table)
#install.packages("lme4")
#install.packages("readxl")
#install.packages("openxlsx")
#install.packages("dplyrx")
#install.packages("tidyr")
#install.packages("Matrix")

data <- read.csv("pheno_height.csv")
split_parts <- strsplit(data$morris, "_", fixed = TRUE)
PI_parts <- sapply(split_parts, `[`, 2)
number_parts <- sub("PI([0-9]+).*", "\\1", PI_parts)
number_parts = paste0("PI", number_parts)
data$PI <- PI_parts
data$number <- number_parts

data_filt = data[, c("number", "DAP_76")]
names(data_filt)[names(data_filt) == "number"] <- "id"



data$Pedigree <- gsub("PI ", "PI", data$Pedigree)
data$plant_height_rep_1<- as.numeric(data$plant_height_rep_1)
data$plant_height_rep_2<- as.numeric(data$plant_height_rep_2)
data$plant_height_rep_3<- as.numeric(data$plant_height_rep_3)
data <- pivot_longer(data, cols = starts_with("Whole_Dry"), names_to = "replicate", values_to = "DW")

data1 <- filter(data, Treatment == "HN")
data2 <- filter(data, Treatment == "LN")


data$Treatment = factor(data$Treatment)
data$Bench = factor(data$Bench)
data$Row = factor(data$Row)
data$Column = factor(data$Column)
data$Pot_id = factor(data$Pot_id)
data$Pedigree = factor(data$Pedigree)
#mod = lm(height ~ Treatment + Pedigree + Treatment*Pedigree, data = data)
#mod = lm(height ~ Pedigree + Row + Column + Pot_id, data = data1)

mod = lmer(DW ~ Pedigree + (1|Row) + (1|Column) + (1|Pot_id), data = data2)
x = summary(mod)
anova = anova(mod)
BLUE <- as.data.frame(fixef(mod))
rownames(BLUE)=gsub("Pedigree","",rownames(BLUE))
BLUE$Pedigree = rownames(BLUE)
BLUE$Pedigree[1] = "PI152651"
HN=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
HN$id[1] = "PI152651"
#write.table(b, file = "data_frame_b.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

LN=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
LN$id[1] = "PI152651"
names(HN)[names(HN) == "Phe"] = "HN"
names(LN)[names(LN) == "Phe"] = "LN"
#write.table(b, file = "data_frame_b.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#map <- read.table("mvp.geno.map" , head = TRUE)
Phenotype = fread("snps.fam")

phenotype1=merge(Phenotype, data_filt, by.x="V1", by.y="id", all.x=T)
num_na_values <- sum(is.na(HN$Phe))
phenotype2=merge(phenotype1, LN, by.x="V1", by.y="id", all.x=T)
phenotypee=phenotype1[,c(1,7)]
names(phenotypee)[names(phenotypee) == "V1"] = "Taxa"
write.table(phenotypee, file = "pheno_height.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
f = read.table("DW.txt")


