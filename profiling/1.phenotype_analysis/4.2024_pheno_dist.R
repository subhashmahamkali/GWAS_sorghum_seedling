data = read.csv("data/1.Phenotype_data/3.2024/1.2024_SAP_Pheno_raw.csv")                            
data = data[,c(1,2,6,7,8,9,11,12,13)]
data$geno <- gsub("PI ", "PI", data$geno)
#node_number = data[,c(1,2,6,10)]
#str(node_number)
ch = data[,c(1,2,3,4,5,6)]
ph = data[,c(1,2,3,7,8,9)]
ph$Plant.Height.1 = ph$Plant.Height.1*100
ph$Plant.Height.2 = ph$Plant.Height.2*100
ph$Plant.Height.3 = ph$Plant.Height.3*100
str(ph)

library(tidyr)
library(dplyr)
library(ggplot2)

ch_wide <- ch %>%
  mutate(Chl_avg = rowMeans(select(., starts_with("Chl")), na.rm = TRUE)) %>%
  group_by(geno, experiment) %>%
  summarise(Chl_avg = mean(Chl_avg, na.rm = TRUE)) %>%
  pivot_wider(names_from = experiment, values_from = Chl_avg)

ch_wide$HN <- as.numeric(ch_wide$HN)
ch_wide$LN <- as.numeric(ch_wide$LN)

# Calculate the means
mean_HN <- mean(ch_wide$HN, na.rm = TRUE)
mean_LN <- mean(ch_wide$LN, na.rm = TRUE)

# Perform paired t-test
t_test <- t.test(ch_wide$HN, ch_wide$LN, paired = TRUE)

# Create density plot
dw <- ggplot(ch_wide) +
  geom_density(aes(x = HN, fill = "HN"), alpha = 0.5) +
  geom_density(aes(x = LN, fill = "LN"), alpha = 0.5) +
  geom_vline(xintercept = mean_HN, linetype = "dashed", color = "red") +
  geom_vline(xintercept = mean_LN, linetype = "dashed", color = "green") +
  scale_fill_manual(values = c("HN" = "red", "LN" = "green")) +
  theme_minimal() +
  labs(
    title = "chlorophyll raw_data 2024",
    x = "chl",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_blank()
  )

# Add significance label
sig_label <- ifelse(
  t_test$p.value < 0.001, '***',
  ifelse(t_test$p.value < 0.01, '**',
         ifelse(t_test$p.value < 0.05, '*', 'ns')
  )
)

x_position <- mean(c(mean_HN, mean_LN))
y_position <- max(c(density(ch_wide$HN, na.rm = TRUE)$y, density(ch_wide$LN, na.rm = TRUE)$y))

# Add the significance label to the plot
dw <- dw +
  annotate("text", x = x_position, y = y_position, label = sig_label, size = 6, color = "blue")

# Display the plot
print(dw)
ggsave(filename = "graphs/1.phenotype_distribution/2024/1.chlorophyll_2024.png", plot = dw,width =5 ,height =3 )




# Step 1: Process the data for Plant Height
ph_wide <- ph %>%
  mutate(Height_avg = rowMeans(select(., starts_with("Plant.Height")), na.rm = TRUE)) %>%
  group_by(geno, experiment) %>%
  summarise(Height_avg = mean(Height_avg, na.rm = TRUE)) %>%
  pivot_wider(names_from = experiment, values_from = Height_avg)

# Convert HN and LN to numeric (if needed)
ph_wide$HN <- as.numeric(ph_wide$HN)
ph_wide$LN <- as.numeric(ph_wide$LN)

# Step 2: Calculate the means
mean_HN <- mean(ph_wide$HN, na.rm = TRUE)
mean_LN <- mean(ph_wide$LN, na.rm = TRUE)

# Step 3: Perform paired t-test
t_test <- t.test(ph_wide$HN, ph_wide$LN, paired = TRUE)

# Step 4: Create the density plot
ph <- ggplot(ph_wide) +
  geom_density(aes(x = HN, fill = "HN"), alpha = 0.5) +
  geom_density(aes(x = LN, fill = "LN"), alpha = 0.5) +
  geom_vline(xintercept = mean_HN, linetype = "dashed", color = "red") +
  geom_vline(xintercept = mean_LN, linetype = "dashed", color = "green") +
  scale_fill_manual(values = c("HN" = "red", "LN" = "green")) +
  theme_minimal() +
  labs(
    title = "Plant Height Raw Data 2024",
    x = "Plant Height (cm)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_blank()
  )

# Step 5: Add significance label
sig_label <- ifelse(
  t_test$p.value < 0.001, '***',
  ifelse(t_test$p.value < 0.01, '**',
         ifelse(t_test$p.value < 0.05, '*', 'ns')
  )
)

x_position <- mean(c(mean_HN, mean_LN))
y_position <- max(c(density(ph_wide$HN, na.rm = TRUE)$y, density(ph_wide$LN, na.rm = TRUE)$y))

ph <- ph +
  annotate("text", x = x_position, y = y_position, label = sig_label, size = 6, color = "blue")

# Step 6: Display the plot
print(ph)

ggsave(filename = "graphs/1.phenotype_distribution/2024/1.plant_height_2024.png", plot = ph,width =5 ,height =3 )









node_number_wide <- node_number %>%
  group_by(geno, experiment) %>%
  summarise(Node_avg = mean(Node.Number, na.rm = TRUE)) %>%
  pivot_wider(names_from = experiment, values_from = Node_avg)

