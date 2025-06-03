library(data.table)
library(ggplot2)

d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_50pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_10pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_1pct.txt")
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_0.1pct.txt")


dim(d)
vcf_header <- readLines("/Users/subhashmahamkali/Downloads/vcf_header.txt")
vcf_individuals <- strsplit(vcf_header[grep("^#CHROM", vcf_header)], "\t")[[1]][-c(1:9)]

# Make sure to assign to gene_load_summary columns (after "chr", "pos", "start", etc.)
colnames(d)[-c(1:4)] <- vcf_individuals
imp = fread("/Users/subhashmahamkali/Downloads/improved.txt")
lan = fread("/Users/subhashmahamkali/Downloads/landrace.txt")
wil = fread("/Users/subhashmahamkali/Downloads/wild.txt")


#imp = imp
#lan = land
#wil = wild

individual_names = imp$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
imp_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns = imp_d[, -c(1:4), with = FALSE]
average_genetic_load = rowMeans(genetic_load_columns, na.rm = TRUE)
imp_d$average_genetic_load = average_genetic_load
average_genetic_load[1]
imp_d$average_genetic_load[1]


individual_names = lan$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
lan_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns = lan_d[, -c(1:4), with = FALSE]
average_genetic_load = rowMeans(genetic_load_columns, na.rm = TRUE)
lan_d$average_genetic_load = average_genetic_load

individual_names = wil$V1
matching_columns = colnames(d)[colnames(d) %in% individual_names]
wil_d = d[, c("geneID", "chr", "start", "end", matching_columns), with = FALSE]
genetic_load_columns = wil_d[, -c(1:4), with = FALSE]
average_genetic_load = rowMeans(genetic_load_columns, na.rm = TRUE)
wil_d$average_genetic_load = average_genetic_load


imp_d$sub_population <- "Improved"
lan_d$sub_population <- "Landrace"
wil_d$sub_population <- "Wild"

combined_data <- rbind(imp_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)],
                       lan_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)],
                       wil_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)])

ggplot(combined_data, aes(x = sub_population, y = average_genetic_load, fill = sub_population)) +
  geom_violin(trim = FALSE) +  # trim = FALSE ensures the violin extends to the full range of data
  theme_minimal() +
  labs(
    title = "Genetic Load Distribution Across Sub-populations",
    x = "Sub-population",
    y = "Average Genetic Load"
  ) +
  scale_fill_manual(values = c("blue", "green", "red")) +  # Customize colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center title
    plot.subtitle = element_text(hjust = 0.5),  # Center subtitle
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.25)
  ) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black")  # Add mean points




landrace_data <- combined_data[sub_population == "Landrace"]
mean(landrace_data$average_genetic_load, na.rm = TRUE)

Improved <- combined_data[sub_population == "Improved"]
mean(Improved$average_genetic_load, na.rm = TRUE)

wild <- combined_data[sub_population == "Wild"]
mean(wild$average_genetic_load, na.rm = TRUE)
library(ggplot2)
ggplot(combined_data, aes(x = sub_population, y = average_genetic_load, fill = sub_population)) +
  geom_boxplot() +
  theme_minimal()
imp_d$num_non_missing <- rowSums(!is.na(genetic_load_columns))
lan_d$num_non_missing <- rowSums(!is.na(genetic_load_columns))
wil_d$num_non_missing <- rowSums(!is.na(genetic_load_columns))
boxplot(imp_d$num_non_missing, lan_d$num_non_missing, wil_d$num_non_missing,
        names = c("Improved", "Landrace", "Wild"),
        ylab = "Non-missing genotype count per gene",
        main = "Coverage per Gene by Group")


indivs <- imp$V1
group_size <- 107
cols <- colnames(d)[colnames(d) %in% indivs]
imp_d <- d[, c("geneID", "chr", "start", "end", cols), with = FALSE]
genetic_load <- imp_d[, -c(1:4), with = FALSE]
imp_d$average_genetic_load <- rowSums(genetic_load, na.rm = TRUE) / group_size
imp_d$sub_population <- "Improved"

# 🟩 Landrace
indivs <- lan$V1
group_size <- 129
cols <- colnames(d)[colnames(d) %in% indivs]
lan_d <- d[, c("geneID", "chr", "start", "end", cols), with = FALSE]
genetic_load <- lan_d[, -c(1:4), with = FALSE]
lan_d$average_genetic_load <- rowSums(genetic_load, na.rm = TRUE) / group_size
lan_d$sub_population <- "Landrace"

# 🟥 Wild
indivs <- wil$V1
group_size <- 50
cols <- colnames(d)[colnames(d) %in% indivs]
wil_d <- d[, c("geneID", "chr", "start", "end", cols), with = FALSE]
genetic_load <- wil_d[, -c(1:4), with = FALSE]
wil_d$average_genetic_load <- rowSums(genetic_load, na.rm = TRUE) / group_size
wil_d$sub_population <- "Wild"

#----------------------------
# Step 5: Combine all data
#----------------------------
combined_data <- rbind(
  imp_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)],
  lan_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)],
  wil_d[, .(geneID, chr, start, end, average_genetic_load, sub_population)]
)

#----------------------------
# Step 6: Violin plot
#----------------------------
five = ggplot(combined_data, aes(x = sub_population, y = average_genetic_load, fill = sub_population)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black") +
  theme_minimal() +
  labs(
    title = "Genetic Load Distribution Across Sub-populations",
    x = "Sub-population",
    y = "Average Genetic Load"
  ) +
  scale_fill_manual(values = c("blue", "green", "red")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.25)
  )
combined_data[, .(mean_load = mean(average_genetic_load)), by = sub_population]


pointone = ggplot(combined_data, aes(x = sub_population, y = average_genetic_load, fill = sub_population)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black") +
  labs(
    title = "Genetic Load Distribution Across Sub-populations top 0.1%",
    x = "Sub-population",
    y = "Average Genetic Load"
  ) +
  theme_minimal()
ggsave("genetic_load_0_1.pdf", width = 7, height = 4)


library(patchwork)
combined_plot <- (five | ten) / (one | pointone)
ggsave("genetic_load_grid.pdf", plot = combined_plot, width = 12, height = 8)
