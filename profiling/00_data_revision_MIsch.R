# Load necessary libraries
library(readxl) # for reading Excel files
library(openxlsx) # for writing Excel files
library(dplyr)
# Step 1: Read the data
csv_data <- read.csv("data/Revised_exp_data/PH_hong.csv") 
excel_data <- read_excel("/Users/subhashmahamkali/Desktop/FW_RE_1.xlsx") 

excel_data_cleaned <- excel_data[ , -c(2, 3, 6, 7, 8, 9)]

csv_data_subset <- csv_data[, c("Pot_id", "Pedigree", "Bench", "Row", "Column", "Treatment")]

# Step 2: Merge the data based on 'Pot_id' - this time focusing on adding the CSV columns to the Excel data
# Assuming 'excel_data' does not initially contain 'ColumnOfInterest', it will be added from 'csv_data'
merged_data <- merge(excel_data_cleaned, csv_data_subset, by = "Pot_id", all.x = TRUE)

# Step 3: Write the merged data to a new Excel file
write.csv(merged_data, file = "data/Revised_exp_data/FW_Revised.csv") # Update the save path
