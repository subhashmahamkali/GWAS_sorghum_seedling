# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)

# Read the Excel file
# Replace 'path_to_your_excel_file.xlsx' with the actual path to your Excel file
df <- read_excel('Copy of PH.xlsx')
df$Pedigree <- gsub("PI ", "PI", df$Pedigree)
# Calculate the mean of the three replicates

df <- df %>%
  mutate(across(starts_with("plant_height_rep_"), ~ as.numeric(.), .names = "{.col}_numeric"))

df <- df %>%
  mutate(MeanDryWeight = rowMeans(select(., ends_with("_numeric")), na.rm = TRUE))

df <- df %>%
  mutate(MeanDryWeight = rowMeans(select(., starts_with("leaf_Count"))))

# View the first few rows of the dataframe to check the MeanDryWeight column
head(df)

# Filter data for High Nitrogen (HN) treatment
df_HN <- df %>% 
  filter(Treatment == "HN")

# Filter data for Low Nitrogen (LN) treatment
df_LN <- df %>% 
  filter(Treatment == "LN")

# Plotting the histogram for HN treatment
ggplot(df_HN, aes(x = MeanDryWeight)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(x = "Mean plant height", y = "Frequency", title = "Distribution of Mean plant height (HN)") +
  theme_minimal()

# Plotting the histogram for LN treatment
ggplot(df_LN, aes(x = MeanDryWeight)) +
  geom_histogram(binwidth = 0.1, fill = "red", color = "black") +
  labs(x = "Mean plant height", y = "Frequency", title = "Distribution of Mean plant height (LN)") +
  theme_minimal()







# Adjusting font size for the high nutrient (HN) treatment histogram
ggplot(df_HN, aes(x = MeanDryWeight)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(x = " plant height(cm)", y = "Frequency", title = "Distribution of average plant height (HN)") +
  theme_minimal() +
  theme(text = element_text(size = 14),  # Base text size for all text elements
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 12),  # Axis text
        plot.title = element_text(size = 18, hjust = 0.5))  # Plot title

# Adjusting font size for the low nutrient (LN) treatment histogram
ggplot(df_LN, aes(x = MeanDryWeight)) +
  geom_histogram(binwidth = 0.1, fill = "red", color = "black") +
  labs(x = "plant height(cm) ", y = "Frequency", title = "Distribution of average plant height (LN)") +
  theme_minimal() +
  theme(text = element_text(size = 14),  # Base text size for all text elements
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 12),  # Axis text
        plot.title = element_text(size = 18, hjust = 0.5))  # Plot title

