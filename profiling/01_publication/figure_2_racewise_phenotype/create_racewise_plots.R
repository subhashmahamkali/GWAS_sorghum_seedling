#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(readr)
  library(stringr)
})

# ---------------------------
# Paths
# ---------------------------
repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"
script_path <- file.path(repo_root, "profiling/01_publication/figure_2_racewise_phenotype/create_racewise_plots.R")
blue_file <- file.path(repo_root, "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021.csv")
cluster_file <- file.path(repo_root, "data/0.SAP/sorted_cluster_data.txt")
out_root <- file.path(repo_root, "graphs/01_publication/figure_2_racewise_phenotype")
script_copy_dir <- file.path(out_root, "scripts")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(script_copy_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helpers
# ---------------------------
normalize_id <- function(x) {
  toupper(gsub("\\s+", "", as.character(x)))
}

safe_name <- function(x) {
  x %>%
    gsub("[^A-Za-z0-9]+", "_", ., perl = TRUE) %>%
    gsub("^_+|_+$", "", ., perl = TRUE)
}

make_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

save_plot <- function(plot_obj, outfile_base, width = 11, height = 7) {
  ggsave(paste0(outfile_base, ".png"), plot_obj, width = width, height = height, dpi = 400, bg = "white")
  ggsave(paste0(outfile_base, ".pdf"), plot_obj, width = width, height = height, bg = "white")
}

# ---------------------------
# Read data
# ---------------------------
blue <- read_csv(blue_file, show_col_types = FALSE)
cluster <- read.table(cluster_file, header = TRUE, stringsAsFactors = FALSE)

cluster_names <- c(
  "1" = "durra",
  "2" = "kafir",
  "3" = "caudatum",
  "4" = "guinea",
  "5" = "milo/durra-bicolor",
  "6" = "mixed/bicolor"
)

# PCA palette (kept consistent across all plots)
pca_colors <- c(
  "durra" = "#FF7F00",
  "kafir" = "#E69F00",
  "caudatum" = "#984EA3",
  "guinea" = "#56B4E9",
  "milo/durra-bicolor" = "#F781BF",
  "mixed/bicolor" = "#A6A6A6"
)

cluster <- cluster %>%
  transmute(
    genotype = normalize_id(PI_numbers),
    ClusterName = cluster_names[as.character(Cluster)]
  ) %>%
  filter(!is.na(ClusterName))

blue <- blue %>%
  mutate(genotype = normalize_id(genotype))

merged <- blue %>%
  inner_join(cluster, by = "genotype")

if (nrow(merged) == 0) {
  stop("No genotype overlap found between BLUE file and cluster file after ID normalization.")
}

# ---------------------------
# Plot functions
# ---------------------------
make_hnln_ridge <- function(df_trait, year, trait) {
  race_order <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ord = median(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ord)) %>%
    pull(ClusterName)

  df_trait <- df_trait %>%
    mutate(ClusterName = factor(ClusterName, levels = race_order))

  mean_df <- df_trait %>%
    group_by(ClusterName, Condition) %>%
    summarise(ref_value = median(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ymid = as.numeric(ClusterName),
      y = ifelse(Condition == "HN", ymid + 0.18, ymid - 0.18),
      yend = ifelse(Condition == "HN", ymid + 0.42, ymid - 0.42)
    )

  ggplot(df_trait, aes(x = Value, y = ClusterName)) +
    geom_density_ridges(
      aes(fill = Condition, group = interaction(ClusterName, Condition)),
      scale = 1.05,
      rel_min_height = 0.01,
      color = "white",
      linewidth = 0.2,
      alpha = 0.62
    ) +
    geom_segment(
      data = mean_df,
      aes(x = ref_value, xend = ref_value, y = y, yend = yend, color = Condition),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.75
    ) +
    scale_fill_manual(values = c("HN" = "#d95f02", "LN" = "#1b9e77")) +
    scale_color_manual(values = c("HN" = "#d95f02", "LN" = "#1b9e77"), guide = "none") +
    labs(
      title = paste0(year, " | ", trait, " | Ridge (HN vs LN)"),
      x = "BLUE value",
      y = "Race",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )
}

make_nr_ridge <- function(df_trait, year, trait) {
  race_order <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ord = median(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ord)) %>%
    pull(ClusterName)

  df_trait <- df_trait %>%
    mutate(ClusterName = factor(ClusterName, levels = race_order))

  mean_df <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ref_value = median(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ymid = as.numeric(factor(ClusterName, levels = levels(df_trait$ClusterName))),
      y = ymid - 0.34,
      yend = ymid + 0.34
    )

  ggplot(df_trait, aes(x = Value, y = ClusterName, fill = ClusterName)) +
    geom_vline(xintercept = 0, color = "#444444", linewidth = 0.45, linetype = "dashed") +
    geom_density_ridges(alpha = 0.75, scale = 1.05, rel_min_height = 0.01, color = "white", linewidth = 0.2) +
    geom_segment(
      data = mean_df,
      aes(x = ref_value, xend = ref_value, y = y, yend = yend),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.75,
      color = "#222222"
    ) +
    scale_fill_manual(values = pca_colors, drop = FALSE) +
    labs(
      title = paste0(year, " | ", trait, " | Ridge (NR)"),
      x = "Nitrogen response (NR)",
      y = "Race",
      fill = "Race"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )
}

make_hnln_violin <- function(df_trait, year, trait) {
  race_order <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ord = median(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ord)) %>%
    pull(ClusterName)

  df_trait <- df_trait %>%
    mutate(ClusterName = factor(ClusterName, levels = race_order))

  ggplot(df_trait, aes(x = ClusterName, y = Value, fill = Condition)) +
    geom_violin(
      aes(group = interaction(ClusterName, Condition)),
      position = position_dodge(width = 0.85),
      linewidth = 0.2,
      trim = TRUE,
      alpha = 0.62
    ) +
    geom_boxplot(
      aes(group = interaction(ClusterName, Condition)),
      position = position_dodge(width = 0.85),
      width = 0.18,
      outlier.size = 0.4,
      alpha = 0.85
    ) +
    scale_fill_manual(values = c("HN" = "#d95f02", "LN" = "#1b9e77")) +
    labs(
      title = paste0(year, " | ", trait, " | Violin+Box (HN vs LN)"),
      x = "Race",
      y = "BLUE value",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

make_nr_violin <- function(df_trait, year, trait) {
  race_order <- df_trait %>%
    group_by(ClusterName) %>%
    summarise(ord = median(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ord)) %>%
    pull(ClusterName)

  df_trait <- df_trait %>%
    mutate(ClusterName = factor(ClusterName, levels = race_order))

  ggplot(df_trait, aes(x = ClusterName, y = Value, fill = ClusterName)) +
    geom_violin(alpha = 0.60, linewidth = 0.2, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.90) +
    scale_fill_manual(values = pca_colors, drop = FALSE) +
    labs(
      title = paste0(year, " | ", trait, " | Violin+Box (NR)"),
      x = "Race",
      y = "Nitrogen response (NR)",
      fill = "Race"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

# ---------------------------
# Build year-wise long data and plot all
# ---------------------------
years <- c("2020", "2021")

for (yr in years) {
  ln_cols <- names(merged)[grepl(paste0("^", yr, ":LN:"), names(merged))]
  hn_cols <- names(merged)[grepl(paste0("^", yr, ":HN:"), names(merged))]
  nr_cols <- names(merged)[grepl(paste0("^", yr, ":.*:NR$"), names(merged))]

  traits_ln <- sub(paste0("^", yr, ":LN:"), "", ln_cols)
  traits_hn <- sub(paste0("^", yr, ":HN:"), "", hn_cols)
  traits_nr <- sub(paste0("^", yr, ":"), "", sub(":NR$", "", nr_cols))
  common_traits <- intersect(intersect(traits_ln, traits_hn), traits_nr)

  if (length(common_traits) == 0) {
    message("No common HN/LN/NR traits for year ", yr, ". Skipping.")
    next
  }

  hnln_cols <- c(
    "genotype", "ClusterName",
    unlist(lapply(common_traits, function(tr) c(paste0(yr, ":HN:", tr), paste0(yr, ":LN:", tr))))
  )
  nr_keep_cols <- c(
    "genotype", "ClusterName",
    unlist(lapply(common_traits, function(tr) paste0(yr, ":", tr, ":NR")))
  )

  long_hnln <- merged %>%
    select(all_of(hnln_cols)) %>%
    pivot_longer(
      cols = -c(genotype, ClusterName),
      names_to = c("Year", "Condition", "Trait"),
      names_sep = ":",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value))

  long_nr <- merged %>%
    select(all_of(nr_keep_cols)) %>%
    pivot_longer(
      cols = -c(genotype, ClusterName),
      names_to = c("Year", "Trait", "Condition"),
      names_pattern = "^([^:]+):(.+):(NR)$",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value))

  # Output directories
  dir_ridge_hnln <- make_dir(out_root, "ridge_hn_ln", yr)
  dir_ridge_nr <- make_dir(out_root, "ridge_nr", yr)
  dir_violin_hnln <- make_dir(out_root, "violin_hn_ln", yr)
  dir_violin_nr <- make_dir(out_root, "violin_nr", yr)

  for (tr in sort(common_traits)) {
    tr_safe <- safe_name(tr)
    df_hnln <- long_hnln %>% filter(Trait == tr)
    df_nr <- long_nr %>% filter(Trait == tr)

    if (nrow(df_hnln) > 0) {
      p_ridge_hnln <- make_hnln_ridge(df_hnln, yr, tr)
      p_violin_hnln <- make_hnln_violin(df_hnln, yr, tr)
      save_plot(p_ridge_hnln, file.path(dir_ridge_hnln, paste0(yr, "_", tr_safe, "_ridge_HN_LN")))
      save_plot(p_violin_hnln, file.path(dir_violin_hnln, paste0(yr, "_", tr_safe, "_violin_box_HN_LN")))
    }

    if (nrow(df_nr) > 0) {
      p_ridge_nr <- make_nr_ridge(df_nr, yr, tr)
      p_violin_nr <- make_nr_violin(df_nr, yr, tr)
      save_plot(p_ridge_nr, file.path(dir_ridge_nr, paste0(yr, "_", tr_safe, "_ridge_NR")))
      save_plot(p_violin_nr, file.path(dir_violin_nr, paste0(yr, "_", tr_safe, "_violin_box_NR")))
    }
  }

  message("Finished year ", yr, " with ", length(common_traits), " common traits.")
}

# Save a copy of this script into output/scripts
if (file.exists(script_path)) {
  file.copy(script_path, file.path(script_copy_dir, "create_racewise_plots.R"), overwrite = TRUE)
}

message("All outputs saved under: ", out_root)
