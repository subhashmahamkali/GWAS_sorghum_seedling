#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(readr)
  library(patchwork)
})

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"
blue_file <- file.path(repo_root, "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021.csv")
cluster_file <- file.path(repo_root, "data/0.SAP/sorted_cluster_data.txt")
out_dir <- file.path(repo_root, "graphs/01_publication/figure_2_racewise_phenotype/combined")
script_copy_dir <- file.path(out_dir, "scripts")
out_pdf <- file.path(out_dir, "2020_NR_plantHeight_daysToFlower_paniclesPerPlot_3cols_A4h21cm.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_copy_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  toupper(gsub("\\s+", "", as.character(x)))
}

cluster_names <- c(
  "1" = "Durra",
  "2" = "Kafir",
  "3" = "Caudatum",
  "4" = "Guinea",
  "5" = "Milo/Durra-\nBicolor",
  "6" = "Bicolor"
)

pca_colors <- c(
  "Durra" = "#F28E2B",
  "Kafir" = "#DCA73C",
  "Caudatum" = "#A56CC1",
  "Guinea" = "#77BCE8",
  "Milo/Durra-\nBicolor" = "#E88AC6",
  "Bicolor" = "#B9B9B9"
)

trait_labels <- c(
  "plantHeight" = "Plant height",
  "daysToFlower" = "Days to flower",
  "paniclesPerPlot" = "Panicles per plot"
)

target_traits <- c("plantHeight", "daysToFlower", "paniclesPerPlot")
nr_cols <- paste0("2020:", target_traits, ":NR")

blue <- read_csv(blue_file, show_col_types = FALSE) %>%
  mutate(genotype = normalize_id(genotype))

cluster <- read.table(cluster_file, header = TRUE, stringsAsFactors = FALSE) %>%
  transmute(
    genotype = normalize_id(PI_numbers),
    ClusterName = cluster_names[as.character(Cluster)]
  ) %>%
  filter(!is.na(ClusterName))

merged <- inner_join(blue, cluster, by = "genotype")

missing_cols <- setdiff(nr_cols, names(merged))
if (length(missing_cols) > 0) {
  stop("Missing NR columns: ", paste(missing_cols, collapse = ", "))
}

long_nr <- merged %>%
  select(genotype, ClusterName, all_of(nr_cols)) %>%
  pivot_longer(
    cols = -c(genotype, ClusterName),
    names_to = c("Year", "Trait", "Condition"),
    names_pattern = "^([^:]+):(.+):(NR)$",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

make_nr_ridge <- function(df_trait, trait_label) {
  race_order <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ord = median(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ord)) %>%
    pull(ClusterName)

  df_trait <- df_trait %>%
    mutate(ClusterName = factor(ClusterName, levels = race_order))

  median_df <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ref_value = median(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ymid = as.numeric(factor(ClusterName, levels = levels(df_trait$ClusterName))),
      y = ymid - 0.34,
      yend = ymid + 0.34
    )

  ggplot(df_trait, aes(x = Value, y = ClusterName, fill = ClusterName)) +
    geom_vline(xintercept = 0, color = "#222222", linewidth = 0.45, linetype = "dotted") +
    geom_density_ridges(alpha = 0.78, scale = 1.05, rel_min_height = 0.01, color = "white", linewidth = 0.18) +
    geom_segment(
      data = median_df,
      aes(x = ref_value, xend = ref_value, y = y, yend = yend),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.5,
      color = "#222222"
    ) +
    scale_fill_manual(values = pca_colors, drop = FALSE) +
    labs(title = paste0(trait_label, " (2020)"), x = "Nitrogen response (NR)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6)
    )
}

plot_list <- lapply(target_traits, function(tr) {
  df <- long_nr %>% filter(Year == "2020", Trait == tr)
  make_nr_ridge(df, trait_labels[[tr]])
})

combined_plot <- wrap_plots(plot_list, ncol = 3, byrow = TRUE)

# Match Figure 1 (A4 landscape) height: 21 cm.
ggsave(out_pdf, combined_plot, width = 29.7, height = 21, units = "cm", bg = "white")

file.copy("/tmp/generate_2020_nr_three_traits_one_row.R",
          file.path(script_copy_dir, "generate_2020_nr_three_traits_one_row.R"),
          overwrite = TRUE)

message("Saved: ", out_pdf)
