library(dplyr)
library(ggplot2)

null_medians_wild_landrace <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/wt_la_null_medians.rds")
merged_median_wild_landrace <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/wt_la_obs_median_final.rds")

null_medians_landrace_improved <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/la_im_null_medians.rds")
merged_median_landrace_improved <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/la_im_obs_median_final.rds")

null_medians_balancing <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/bnull_medians.rds")
merged_median_balancing <- readRDS("/Users/subhashmahamkali/Downloads/gwas_sap/bmerged_median_final.rds")

merged_median_wild_landrace
#0.09396357
merged_median_landrace_improved
#0.09490959
merged_median_balancing
#0.1071778


data_combined <- data.frame(
  group = rep(c("Wild vs Landrace", "Landrace vs Improved", "Balancing Selection"), 
              c(length(null_medians_wild_landrace), 
                length(null_medians_landrace_improved), 
                length(null_medians_balancing))),
  zs = c(null_medians_wild_landrace, null_medians_landrace_improved, null_medians_balancing)
)

x_positions <- c("Balancing Selection" = 1, "Landrace vs Improved" = 2, "Wild vs Landrace" = 3)

p_combined <- ggplot(data_combined, aes(x = group, y = zs)) +
  geom_violin(aes(fill = group), color = "black", alpha = 0.3) +
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.7, size = 2) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = merged_median_balancing, yend = merged_median_balancing),
               color = "red", linetype = "dashed", size = 1) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = merged_median_landrace_improved, yend = merged_median_landrace_improved),
               color = "red", linetype = "dashed", size = 1) +
  geom_segment(aes(x = 2.8, xend = 3.2, y = merged_median_wild_landrace, yend = merged_median_wild_landrace),
               color = "red", linetype = "dashed", size = 1) +
  labs(title = "",
       y = "Median Zero-Shot Score", x = "") +
  annotate("text", x = 3, y = merged_median_wild_landrace, 
           label = paste("Observed value = ",round(merged_median_wild_landrace, 3)),
           color = "red", vjust = -1, fontface = "bold") +
  annotate("text", x = 2, y = merged_median_landrace_improved, 
           label = paste("Observed value = ",round(merged_median_landrace_improved, 3)),
           color = "red", vjust = -1, fontface = "bold") +
  annotate("text", x = 1, y = merged_median_balancing, 
           label = paste("Observed value = ",round(merged_median_balancing, 3)),
           color = "red", vjust = -1, fontface = "bold") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),   
    axis.text.x  = element_text(size = 14),                  
    axis.text.y  = element_text(size = 14)                   
  )
ggsave("/Users/subhashmahamkali/Documents/gwas_sap/graphs/violin_plot.png", p_combined, width = 8, height = 6, dpi = 300)
