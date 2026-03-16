#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(readr)
})

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"
blue_file <- file.path(repo_root, "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021.csv")
cluster_file <- file.path(repo_root, "data/0.SAP/sorted_cluster_data.txt")
out_dir <- file.path(repo_root, "graphs/01_publication/figure_2_racewise_phenotype/individual_nr_a4height")
script_out <- file.path(out_dir, "export_individual_nr_a4height.R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  "branchInternodeLength" = "Branch internode length",
  "daysToFlower" = "Days to flower",
  "extantLeafNumber" = "Extant leaf number",
  "flagLeafLength" = "Flag leaf length",
  "flagLeafWidth" = "Flag leaf width",
  "leafAngleStandardDeviation" = "Leaf angle standard deviation",
  "medianLeafAngle" = "Median leaf angle",
  "panicleGrainWeight" = "Panicle grain weight",
  "plantHeight" = "Plant height",
  "primaryBranchNumber" = "Primary branch number",
  "rachisDiameterLower" = "Rachis diameter lower",
  "rachisDiameterUpper" = "Rachis diameter upper",
  "rachisLength" = "Rachis length",
  "stemDiameterLower" = "Stem diameter lower",
  "stemDiameterUpper" = "Stem diameter upper",
  "thirdLeafLength" = "Third leaf length",
  "thirdLeafWidth" = "Third leaf width",
  "tillersPerPlant" = "Tillers per plant"
)

all_traits_2020 <- c(
  "daysToFlower", "medianLeafAngle", "leafAngleStandardDeviation", "paniclesPerPlot",
  "panicleGrainWeight", "estimatedPlotYield", "flagLeafLength", "flagLeafWidth",
  "extantLeafNumber", "plantHeight", "thirdLeafLength", "thirdLeafWidth",
  "tillersPerPlant", "stemDiameterLower", "stemDiameterUpper", "rachisLength",
  "rachisDiameterLower", "rachisDiameterUpper", "primaryBranchNumber", "branchInternodeLength"
)

exclude_traits <- c("plantHeight", "daysToFlower", "estimatedPlotYield", "paniclesPerPlot")
keep_traits_2020 <- setdiff(all_traits_2020, exclude_traits)
all_traits_2021 <- c("plantHeight", "daysToFlower", "tillersPerPlant", "stemDiameterLower", "stemDiameterUpper")
keep_traits_2021 <- setdiff(all_traits_2021, c("tillersPerPlant"))

blue <- read_csv(blue_file, show_col_types = FALSE) %>%
  mutate(genotype = normalize_id(genotype))

cluster <- read.table(cluster_file, header = TRUE, stringsAsFactors = FALSE) %>%
  transmute(
    genotype = normalize_id(PI_numbers),
    ClusterName = cluster_names[as.character(Cluster)]
  ) %>%
  filter(!is.na(ClusterName))

merged <- inner_join(blue, cluster, by = "genotype")

nr_cols_2020 <- paste0("2020:", keep_traits_2020, ":NR")
nr_cols_2021 <- paste0("2021:", keep_traits_2021, ":NR")

long_nr <- merged %>%
  select(genotype, ClusterName, all_of(c(nr_cols_2020, nr_cols_2021))) %>%
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
    labs(title = trait_label, x = "Nitrogen response (NR)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 12, 12, 12)
    )
}

specs <- tibble::tibble(
  Year = c(rep("2020", length(keep_traits_2020)), rep("2021", length(keep_traits_2021))),
  Trait = c(keep_traits_2020, keep_traits_2021)
) %>%
  mutate(TraitLabel = unname(trait_labels[Trait])) %>%
  filter(!is.na(TraitLabel)) %>%
  arrange(Year, TraitLabel)

# Match Figure 1 A4 landscape height
plot_height_cm <- 21
plot_width_cm <- 19

for (i in seq_len(nrow(specs))) {
  yr <- specs$Year[[i]]
  tr <- specs$Trait[[i]]
  tl <- specs$TraitLabel[[i]]
  label <- paste0(tl, " (", yr, ")")
  df <- long_nr %>% filter(Year == yr, Trait == tr)
  if (nrow(df) == 0) next

  p <- make_nr_ridge(df, label)
  outfile <- file.path(out_dir, paste0(yr, "_", tr, "_NR_ridge_A4h21cm.pdf"))
  ggsave(outfile, p, width = plot_width_cm, height = plot_height_cm, units = "cm", bg = "white")
}

file.copy("/tmp/export_individual_nr_a4height.R", script_out, overwrite = TRUE)
cat("Saved NR individual plots to:", out_dir, "\n")
