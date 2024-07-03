library(readxl)
data = read_excel("SAPMerged2021_v2.3.xlsx")
data2 <- data[, c(2, 7, 11)]
#data2 = na.omit(data2)
df_HN <- data2 %>% 
  filter(Treatment == "HN")
df_LN <- data2 %>% 
  filter(Treatment == "LN")

data_HN <- aggregate(PlantHeight ~ PINumber, data = df_HN, FUN = mean)
data_LN <- aggregate(PlantHeight ~ PINumber, data = df_LN, FUN = mean)

df_HN <- data2 %>% filter(Treatment == "HN")
df_LN <- data2 %>% filter(Treatment == "LN")

names(data_HN)[names(data_HN) == "PlantHeight"] <- "PlantHeight_HN"
names(data_LN)[names(data_LN) == "PlantHeight"] <- "PlantHeight_LN"

data_HN <- data_HN %>%
  rename_with(~ paste0(., "_HN"), -PINumber)
data_LN <- data_LN %>%
  rename_with(~ paste0(., "_LN"), -PINumber)

merged_data <- merge(data_HN, data_LN, by = "PINumber", all = TRUE)
merged_data = na.omit(merged_data)
merged_data$NR = (merged_data$PlantHeight_HN - merged_data$PlantHeight_LN) / (merged_data$PlantHeight_LN)

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

