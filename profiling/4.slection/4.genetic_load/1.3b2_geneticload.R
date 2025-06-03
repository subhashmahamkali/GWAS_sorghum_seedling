library(data.table)

b = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2_gene.bed")
b = b[,c(1:4,9)]
colnames(b)[5] <- "geneID"
d = fread("/Users/subhashmahamkali/Downloads/gene_load_summary_top_1pct.txt")

d$geneID_clean <- sub("ID=(Sobic\\.\\d+G\\d+)\\..*", "\\1", d$geneID)
d$geneID <- d$geneID_clean


lan <- fread("/Users/subhashmahamkali/Downloads/landrace.txt")$V1
lan_ids <- as.character(rownames(as.data.frame(lan)))  # make sure they're character type
lan_cols <- paste0("genetic_load_", lan_ids)

lan_cols <- lan_cols[lan_cols %in% colnames(d)]
d$mean_load_lan <- rowMeans(d[, ..lan_cols], na.rm = TRUE)

b2_merged <- merge(b, d[, .(geneID, mean_load_lan)], by = "geneID")

library(ggplot2)

ggplot(b2_merged, aes(x = V4, y = mean_load_lan)) +
  geom_point(alpha = 0.5) +
  labs(title = "B2 Selection Signal vs Genetic Load (Landrace)",
       x = "B2 Score (V4 column)", y = "Mean Genetic Load (Landraces)") +
  theme_minimal()

b2_threshold <- 1500
load_threshold <- quantile(b2_merged$mean_load_lan, 0.99, na.rm = TRUE) 

ggplot(b2_merged, aes(x = V4, y = mean_load_lan)) +
  geom_point(alpha = 0.4, color = "black") +
  geom_vline(xintercept = b2_threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = load_threshold, linetype = "dashed", color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "B2 Selection vs Genetic Load (Landraces)",
    x = "B2 Score",
    y = "Mean Genetic Load (Landraces)"
  ) +
  annotate("text", x = b2_threshold + 50, y = max(b2_merged$mean_load_lan, na.rm = TRUE),
           label = "B2", color = "blue", hjust = 0) +
  annotate("text", x = min(b2_merged$V4, na.rm = TRUE), y = load_threshold + 10,
           label = "Top 1% Load", color = "red", hjust = 0)


high_b2_high_load <- b2_merged[V4 >= b2_threshold & mean_load_lan >= load_threshold]
fwrite(high_b2_high_load, "high_b2_high_load_genes.txt", sep = "\t")

library(readxl)
a=read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/7.sorghum.annotations.xlsx", col_names = F)
colnames(a) <- a[3, ]
a <- a[-(1:3), ]  
a = a[,c(1,9,22)]

a$Gene <- as.character(a$Gene)
high_b2_high_load$geneID <- as.character(high_b2_high_load$geneID)
a_unique <- a[!duplicated(a$Gene), ]

# Merge to get annotations
annotated_genes <- merge(high_b2_high_load,a_unique[, c("Gene", "Description", "PFAMs")],
  by.x = "geneID",
  by.y = "Gene",
  all.x = TRUE)

# Add a column indicating whether each gene is in the top-right quadrant
b2_merged$highlight <- ifelse(
  b2_merged$V4 >= 1500 & b2_merged$mean_load_lan >= quantile(b2_merged$mean_load_lan, 0.99, na.rm = TRUE),
  "High B2 & High Load",
  "Other")

ggplot(b2_merged, aes(x = V4, y = mean_load_lan, color = highlight)) +
  geom_point(alpha = 0.6, size = 2) +
  # Threshold lines
  geom_vline(xintercept = b2_threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = load_threshold, linetype = "dashed", color = "red", size = 1) +
  # Labels for lines
  annotate("text", x = b2_threshold + 50, y = max(b2_merged$mean_load_lan, na.rm = TRUE),
           label = "B2 Threshold = 1500", color = "blue", angle = 90, vjust = -0.5, size = 3.5) +
  annotate("text", x = min(b2_merged$V4, na.rm = TRUE), y = load_threshold + 10,
           label = "Top 1% Genetic Load", color = "red", hjust = 0, size = 3.5) +
  # Styling
  scale_color_manual(values = c("High B2 & High Load" = "orange", "Other" = "grey50")) +
  theme_minimal() +
  labs(
    title = "B2 vs Genetic Load",
    x = "B2 Score",
    y = "Mean Genetic Load (Landraces)",
    color = "Category"
  )
