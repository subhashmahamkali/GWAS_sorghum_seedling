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

trait_categories <- c(
  "branchInternodeLength" = "Architecture",
  "plantHeight" = "Architecture",
  "primaryBranchNumber" = "Architecture",
  "tillersPerPlant" = "Architecture",
  "daysToFlower" = "Developmental",
  "extantLeafNumber" = "Developmental",
  "flagLeafLength" = "Developmental",
  "flagLeafWidth" = "Developmental",
  "leafAngleStandardDeviation" = "Developmental",
  "medianLeafAngle" = "Developmental",
  "stemDiameterLower" = "Developmental",
  "stemDiameterUpper" = "Developmental",
  "thirdLeafLength" = "Developmental",
  "thirdLeafWidth" = "Developmental",
  "panicleGrainWeight" = "Panicle",
  "rachisDiameterLower" = "Panicle",
  "rachisDiameterUpper" = "Panicle",
  "rachisLength" = "Panicle"
)

all_traits_2020 <- c(
  "daysToFlower", "medianLeafAngle", "leafAngleStandardDeviation", "paniclesPerPlot",
  "panicleGrainWeight", "estimatedPlotYield", "flagLeafLength", "flagLeafWidth",
  "extantLeafNumber", "plantHeight", "thirdLeafLength", "thirdLeafWidth",
  "tillersPerPlant", "stemDiameterLower", "stemDiameterUpper", "rachisLength",
  "rachisDiameterLower", "rachisDiameterUpper", "primaryBranchNumber", "branchInternodeLength"
)

exclude_traits <- c("plantHeight", "daysToFlower", "estimatedPlotYield", "paniclesPerPlot")
keep_traits <- setdiff(all_traits_2020, exclude_traits)
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

nr_cols_2020 <- paste0("2020:", keep_traits, ":NR")
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
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4)
    )
}

plot_specs <- tibble::tibble(
  Year = c(rep("2020", length(keep_traits)), rep("2021", length(keep_traits_2021))),
  Trait = c(keep_traits, keep_traits_2021)
) %>%
  mutate(
    Category = unname(trait_categories[Trait]),
    TraitLabel = unname(trait_labels[Trait])
  ) %>%
  filter(!is.na(Category), !is.na(TraitLabel)) %>%
  arrange(Category, Year, TraitLabel)

category_names <- sort(unique(plot_specs$Category))
category_tags <- setNames(LETTERS[seq_along(category_names)], category_names)

category_plots <- lapply(category_names, function(category_name) {
  category_specs <- plot_specs %>% filter(Category == category_name)

  plot_list <- lapply(seq_len(nrow(category_specs)), function(i) {
    yr <- category_specs$Year[[i]]
    tr <- category_specs$Trait[[i]]
    label <- paste0(category_specs$TraitLabel[[i]], " (", yr, ")")
    make_nr_ridge(filter(long_nr, Year == yr, Trait == tr), label)
  })

  section_title <- wrap_elements(
    full = grid::textGrob(
      paste(category_tags[[category_name]], category_name),
      x = 0.02,
      hjust = 0,
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
  )

  section_body <- wrap_plots(plot_list, ncol = 4, byrow = TRUE)

  section_title / section_body + plot_layout(heights = c(0.12, 1))
})

category_heights <- vapply(category_names, function(category_name) {
  ceiling(sum(plot_specs$Category == category_name) / 4)
}, numeric(1))

combined_plot <- wrap_plots(category_plots, ncol = 1, heights = category_heights)

out_pdf <- file.path(out_dir, "2020_2021_NR_remaining_traits_ridge_4cols.pdf")
ggsave(out_pdf, combined_plot, width = 18, height = 16.5, units = "in", bg = "white")

file.copy(
  "/Users/subhashmahamkali/Documents/gwas_sap/profiling/01_publication/figure_2_racewise_phenotype/combine_2020_nr_remaining_traits.R",
  file.path(script_copy_dir, "combine_2020_nr_remaining_traits.R"),
  overwrite = TRUE
)

message("Saved: ", out_pdf)
