#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_file <- "/Users/subhashmahamkali/Documents/gwas_sap/data/6.genetic_load_sap/miss50/3.combined/sap_genetic_load_by_sample_miss50.tsv"
cluster_file <- "/Users/subhashmahamkali/Documents/gwas_sap/data/0.SAP/sorted_cluster_data.txt"
out_dir <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/3.genetic_load/sap_miss50_boxplots_by_cutoff"
script_copy <- file.path(out_dir, "plot_sap_miss50_boxplots_by_cutoff.R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_names <- c(
  "1" = "Durra",
  "2" = "Kafir",
  "3" = "Caudatum",
  "4" = "Guinea",
  "5" = "Milo/Durra-Bicolor",
  "6" = "Bicolor"
)

cutoff_levels <- c("1%", "5%", "10%", "25%")
race_levels <- c("Bicolor", "Caudatum", "Durra", "Guinea", "Kafir", "Milo/Durra-Bicolor")

dt <- fread(in_file)
cl <- fread(cluster_file)

dt[, PI := sub("_.*$", "", sample)]
cl[, `:=`(
  PI = gsub("\\s+", "", PI_numbers),
  Cluster = as.character(Cluster),
  Race = cluster_names[as.character(Cluster)]
)]

plot_dt <- merge(dt, cl[, .(PI, Race)], by = "PI", all.x = TRUE)
plot_dt <- plot_dt[!is.na(Race)]
plot_dt[, cutoff := factor(cutoff, levels = cutoff_levels)]
plot_dt[, Race := factor(Race, levels = race_levels)]

for (cut in cutoff_levels) {
  d <- plot_dt[cutoff == cut]
  if (nrow(d) == 0) next

  p <- ggplot(d, aes(x = Race, y = Genetic_load, fill = Race)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5) +
    labs(
      title = paste0("SAP Genetic Load (miss50) - ", cut, " cutoff"),
      x = NULL,
      y = "Genetic load (more negative = higher burden)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold")
    )

  out_stub <- paste0("sap_genetic_load_boxplot_", gsub("%", "pct", cut), "_miss50")
  ggsave(file.path(out_dir, paste0(out_stub, ".pdf")), p, width = 8, height = 5, units = "in", bg = "white")
  ggsave(file.path(out_dir, paste0(out_stub, ".png")), p, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
}

combined <- ggplot(plot_dt, aes(x = Race, y = Genetic_load, fill = Race)) +
  geom_boxplot(width = 0.7, outlier.size = 0.4) +
  facet_wrap(~ cutoff, nrow = 2, scales = "free_y") +
  labs(
    title = "SAP Genetic Load by Race Across Cutoffs (miss50)",
    x = NULL,
    y = "Genetic load (more negative = higher burden)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "sap_genetic_load_boxplot_all_cutoffs_miss50.pdf"),
       combined, width = 12, height = 8, units = "in", bg = "white")
ggsave(file.path(out_dir, "sap_genetic_load_boxplot_all_cutoffs_miss50.png"),
       combined, width = 12, height = 8, units = "in", dpi = 300, bg = "white")

file.copy("/tmp/plot_sap_miss50_boxplots_by_cutoff.R", script_copy, overwrite = TRUE)
cat("Saved plots to:", out_dir, "\n")
