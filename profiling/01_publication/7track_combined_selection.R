#!/usr/bin/env Rscript
# ============================================================
# Top panel: Sorghum | Maize | Sunflower at 5% cutoff
# Publication-ready violin + boxplot + raw points
# Fixes sunflower violin compression near zero
# ============================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(patchwork)

# ============================================================
# 1. EDIT THESE PATHS
# ============================================================

sorghum_file   <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/3.genetic_load/genome_level_1_5_10_violin_custom_palette/sorghum/sorghum_genome_level_by_sample_1_5_10.tsv"
maize_file     <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/3.genetic_load/genome_level_1_5_10_violin_custom_palette/maize/maize_genome_level_by_sample_1_5_10.tsv"
sunflower_file <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/3.genetic_load/genome_level_1_5_10_violin_custom_palette/sunflower/sunflower_genome_level_by_sample_1_5_10.tsv"
output_dir     <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/3.genetic_load/combined_A4_panels/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. SETTINGS
# ============================================================

group_levels <- c("wild", "landrace", "improved")

group_colors <- c(
  wild     = "#4DAF4A",
  landrace = "#377EB8",
  improved = "#984EA3"
)

# Wilcoxon pairwise comparisons
comparisons <- list(
  c("wild", "landrace"),
  c("wild", "improved"),
  c("landrace", "improved")
)

# Samples to exclude (reference artifacts)
# Maize: B73 triplicates
# Sunflower: V193 = failed QC, only 85K sites
exclude <- list(
  Maize     = c("282set_B73", "B73", "german_B73"),
  Sunflower = c("V193")
)

# ============================================================
# 3. LOAD DATA
# ============================================================

sorghum   <- fread(sorghum_file)
maize     <- fread(maize_file)
sunflower <- fread(sunflower_file)

# ============================================================
# 4. MAKE ONE PANEL
# ============================================================

make_panel <- function(dt, title, cutoff = "5%", exclude_samples = NULL) {
  
  # Filter to chosen cutoff
  d <- dt[cutoff_lab == cutoff]
  
  # Remove known artifacts
  if (!is.null(exclude_samples)) {
    d <- d[!sample %in% exclude_samples]
  }
  
  # Remove missing values just in case
  d <- d[!is.na(Abs_genetic_load) & !is.na(group)]
  
  # Set group order
  d[, group := factor(group, levels = group_levels)]
  
  # Drop unused factor levels
  d <- droplevels(d)
  
  # Sample sizes for x-axis labels
  n_labels <- d[, .N, by = group]
  n_vec <- setNames(n_labels$N, as.character(n_labels$group))
  x_labels <- setNames(
    paste0(group_levels, "\n(n=", ifelse(group_levels %in% names(n_vec), n_vec[group_levels], 0), ")"),
    group_levels
  )
  
  # Dynamic y-axis range
  y_min <- min(d$Abs_genetic_load, na.rm = TRUE)
  y_max <- max(d$Abs_genetic_load, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Padding so violin is not glued to bottom
  # keep lower limit >= 0 because genetic load cannot be negative
  lower_lim <- max(0, y_min - 0.08 * y_range)
  upper_lim <- y_max + 0.18 * y_range
  
  # Bracket positions
  bracket_base <- y_max + 0.04 * y_range
  bracket_step <- 0.07 * y_range
  
  ggplot(d, aes(x = group, y = Abs_genetic_load, fill = group)) +
    
    # Violin
    geom_violin(
      trim      = TRUE,
      width     = 0.85,
      alpha     = 0.90,
      color     = "grey30",
      linewidth = 0.3
    ) +
    
    # Boxplot
    geom_boxplot(
      width         = 0.14,
      fill          = "white",
      color         = "grey20",
      linewidth     = 0.35,
      outlier.shape = NA
    ) +
    
    # Raw points for transparency
    geom_jitter(
      width  = 0.08,
      size   = 0.9,
      alpha  = 0.35,
      color  = "black",
      stroke = 0
    ) +
    
    # Wilcoxon significance
    stat_compare_means(
      comparisons  = comparisons,
      method       = "wilcox.test",
      label        = "p.signif",
      hide.ns      = TRUE,
      size         = 3.2,
      tip.length   = 0.01,
      bracket.size = 0.4,
      label.y      = c(
        bracket_base,
        bracket_base + bracket_step,
        bracket_base + 2 * bracket_step
      )
    ) +
    
    scale_fill_manual(values = group_colors) +
    
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    
    scale_y_continuous(
      labels = label_comma(),
      limits = c(lower_lim, upper_lim),
      expand = expansion(mult = c(0, 0))
    ) +
    
    labs(
      title = title,
      y     = "Absolute genetic load",
      x     = NULL
    ) +
    
    coord_cartesian(clip = "off") +
    
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.title.y       = element_text(size = 11, face = "bold"),
      axis.text.x        = element_text(size = 10, face = "bold", color = "grey20"),
      axis.text.y        = element_text(size = 9, color = "grey30"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey60", linewidth = 0.5),
      plot.margin        = margin(8, 10, 8, 10)
    )
}

# ============================================================
# 5. BUILD THREE PANELS
# ============================================================

p1 <- make_panel(sorghum,   "Sorghum 5%")
p2 <- make_panel(maize,     "Maize 5%",     exclude_samples = exclude$Maize)
p3 <- make_panel(sunflower, "Sunflower 5%", exclude_samples = exclude$Sunflower)

# Combine side by side
top_panel <- p1 + p2 + p3 +
  plot_annotation(
    title   = "Genome-Wide Absolute Genetic Load (Top 5% Deleterious SNPs)",
    caption = "Wilcoxon rank-sum test; **** p < 0.0001, ** p < 0.01; ns hidden",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5,
                                  margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, color = "grey40", hjust = 0)
    )
  )

# ============================================================
# 6. SAVE
# ============================================================

ggsave(
  file.path(output_dir, "top_panel_5pct_publication.pdf"),
  top_panel,
  width  = 10.5,
  height = 4.8,
  units  = "in",
  bg     = "white",
  useDingbats = FALSE
)

ggsave(
  file.path(output_dir, "top_panel_5pct_publication.png"),
  top_panel,
  width  = 10.5,
  height = 4.8,
  units  = "in",
  dpi    = 600,
  bg     = "white"
)

message("Done! Files saved to: ", output_dir)