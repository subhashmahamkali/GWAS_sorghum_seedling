library(readxl)
library(dplyr)
library(data.table)
data = read_excel("SAP2020_merged_v2.xls")
data2 <- data[, c(3, 10, 31)]
df_HN <- data2 %>% 
  filter(Treatment == "SufficientNitrogen")
df_LN <- data2 %>% 
  filter(Treatment == "LowNitrogen")

data_HN <- aggregate(BranchInternodeLength ~ SorghumAccessionYellow, data = df_HN, FUN = mean)
data_LN <- aggregate(BranchInternodeLength ~ SorghumAccessionYellow, data = df_LN, FUN = mean)

names(data_HN)[names(data_HN) == "BranchInternodeLength"] <- "BranchInternodeLength_HN"
names(data_LN)[names(data_LN) == "BranchInternodeLength"] <- "BranchInternodeLength_LN"

merged_data <- merge(data_HN, data_LN, by = "SorghumAccessionYellow", all = TRUE)
merged_data = na.omit(merged_data)
merged_data$NR = (merged_data$BranchInternodeLength_HN - merged_data$BranchInternodeLength_LN) / (merged_data$BranchInternodeLength_LN)
names(merged_data)[names(merged_data) == "NR"] <- "NR_BranchInternodeLength"
write.table(merged_data, file = "BranchInternodeLength.txt", sep = "\t", row.names = FALSE, col.names = TRUE)



plot(density(merged_data$NR), 
     main = "Distribution of NR", 
     xlab = "NR", 
     ylab = "Density", 
     col = "blue")
min_NR <- round(min(merged_data$NR), 2)
max_NR <- round(max(merged_data$NR), 2)

# Create a density plot of the NR column
p = ggplot(merged_data, aes(x = NR)) +
  geom_density(fill = "gold", alpha = 0.5) +
  labs(title = "2021_Panicles Per Plant", x = "Nitrogen Response", y = "Density") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(min_NR, max_NR, by = 2)) +
  theme_minimal() + theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )
ggsave("PaniclesPerPlant_2021_NR.png", plot = p, width = 8, height = 6)
p
merged_data$TillersPerPlant_HN = merged_data$TillersPerPlant_HN + 0.1
merged_data$TillersPerPlant_LN = merged_data$TillersPerPlant_LN + 0.1


data1 = read.table("BranchInternodeLength.txt", header = T)
data2 = read.table("PrimaryBranchNo.txt", header = T)
data3 = read.table("RachisDiameterUpper.txt", header = T)
data4 = read.table("RachisDiameterLower.txt", header = T)
data5 = read.table("RachisLength.txt", header = T)
data6 = read.table("StemDiameterUpper.txt", header = T)
data7 = read.table("StemDiameterLower.txt", header = T)
data8 = read.table("TillersPerPlant.txt", header = T)
data9 = read.table("ThirdLeafWidth.txt", header = T)
data10 = read.table("ThirdLeafLength.txt", header = T)
data11 = read.table("PlantHeight.txt", header = T)
data12 = read.table("ExtantLeafNumber.txt", header = T)
data13 = read.table("FlagLeafWidth.txt", header = T)
data14 = read.table("FlagLeafLength.txt", header = T)
data15 = read.table("PanicleGrainWeight.txt", header = T)
data16 = read.table("PaniclesPerPlot.txt", header = T)
data17 = read.table("MedianLeafAngle.txt", header = T)
data18 = read.table("DaysToBloom.txt", header = T)
data_list = list(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15,data16, data17, data18)
merged_data = Reduce(function(x, y) merge(x, y, by = "SorghumAccessionYellow", all = TRUE), data_list)
merged_data$SorghumAccessionYellow = gsub(" ", "", merged_data$SorghumAccessionYellow)
names(merged_data)[names(merged_data) == "SorghumAccessionYellow"] = "Taxa"

Phenotype = fread("snps.fam")
GWAS=merge(Phenotype, merged_data, by.x="V1", by.y="Taxa", all.x=T)
GWAS1=GWAS[,c(1,7:60)]
names(GWAS1)[names(GWAS1) == "V1"] = "Taxa"
write.table(GWAS1, file = "2020.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
f = read.table("2020.txt")

