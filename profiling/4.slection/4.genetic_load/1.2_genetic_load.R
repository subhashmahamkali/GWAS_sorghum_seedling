library(data.table)
library(ggplot2)

# Load gene-level genetic load file (e.g., for top 10% deleterious SNPs)
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_50pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_10pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_1pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_0.1pct.txt")


# Load VCF header to get individual IDs
vcf_header <- readLines("/Users/subhashmahamkali/Downloads/vcf_header.txt")
vcf_individuals <- strsplit(vcf_header[grep("^#CHROM", vcf_header)], "\t")[[1]][-c(1:9)]
colnames(d)[-c(1:4)] <- vcf_individuals  # Assign correct individual IDs to genetic load columns

# Load population lists
imp <- fread("/Users/subhashmahamkali/Downloads/improved.txt")
lan <- fread("/Users/subhashmahamkali/Downloads/landrace.txt")
wil <- fread("/Users/subhashmahamkali/Downloads/wild.txt")

### Function to compute average load using full denominator
process_population <- function(pop_df, pop_name, d) {
  indiv_ids <- pop_df$V1
  match_cols <- colnames(d)[colnames(d) %in% indiv_ids]
  
  df_sub <- d[, c("geneID", "chr", "start", "end", match_cols), with = FALSE]
  load_matrix <- df_sub[, -c(1:4), with = FALSE]
  
  # Sum across individuals and divide by total individuals (to penalize missing)
  df_sub$genetic_load_sum <- rowSums(load_matrix, na.rm = TRUE)
  df_sub$average_genetic_load <- df_sub$genetic_load_sum / length(indiv_ids)
  
  df_sub$sub_population <- pop_name
  df_sub$num_non_missing <- rowSums(!is.na(load_matrix))
  
  return(df_sub[, .(geneID, chr, start, end, average_genetic_load, num_non_missing, sub_population)])
}

imp_d <- process_population(imp, "Improved", d)
lan_d <- process_population(lan, "Landrace", d)
wil_d <- process_population(wil, "Wild", d)

# Combine for plotting
combined_data <- rbind(imp_d, lan_d, wil_d)

# OPTIONAL: Filter genes with poor coverage in wild
# combined_data <- combined_data[!(sub_population == "Wild" & num_non_missing < 40)]

# Violin plot of average genetic load
ggplot(combined_data, aes(x = sub_population, y = average_genetic_load, fill = sub_population)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black") +
  theme_minimal() +
  labs(
    title = "Normalized Genetic Load Across Subpopulations",
    x = "Sub-population",
    y = "Normalized Average Genetic Load"
  ) +
  scale_fill_manual(values = c("blue", "green", "red")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Print mean values for interpretation
cat("Mean genetic load by group:\n")
print(combined_data[, .(mean_load = mean(average_genetic_load, na.rm = TRUE)), by = sub_population])

# OPTIONAL: Compare coverage across groups
boxplot(average_genetic_load ~ sub_population, data = combined_data,
        main = "Normalized Genetic Load by Group", ylab = "Avg Load (Normalized)")



individual_names = imp$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
group_size = length(matching_columns)  # Total individuals in group
imp_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns_imp = imp_d[, -c(1:4), with = FALSE]
imp_d$average_genetic_load = rowSums(genetic_load_columns, na.rm = TRUE) / group_size
sum(is.na(genetic_load_columns_imp))


individual_names = lan$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
group_size = length(matching_columns)
lan_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns = lan_d[, -c(1:4), with = FALSE]
lan_d$average_genetic_load = rowSums(genetic_load_columns, na.rm = TRUE) / group_size

individual_names = wil$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
group_size = length(matching_columns)
wil_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns = wil_d[, -c(1:4), with = FALSE]
wil_d$average_genetic_load = rowSums(genetic_load_columns, na.rm = TRUE) / group_size





