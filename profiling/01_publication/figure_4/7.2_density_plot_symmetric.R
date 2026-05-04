library(data.table)
library(ggplot2)
library(patchwork)

#-----------------------------
# paths
#-----------------------------
BASE_DIR_B2  <- "data/4.GWAS_selection/A.balancing_selection"
BASE_DIR_POS <- "data/4.GWAS_selection/B.positive_selection"
OUTDIR       <- "graphs/01_publication/4.GWAS_selection_integration/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# only 3 trait categories
#-----------------------------
traits <- c("Architecture", "Developmental", "Panicle")

#-----------------------------
# original colors
#-----------------------------
cols_b2 <- c(
  Wild     = "#1b9e77",
  Landrace = "#d95f02",
  Improved = "#7570b3"
)

cols_pos <- c(
  Landrace = "#E69F00",
  Improved = "#009E73"
)

#-----------------------------
# helper function: B2 plot
#-----------------------------
make_b2_plot <- function(trait) {
  
  file <- file.path(
    BASE_DIR_B2,
    paste0("NR_", trait, "_B2_top1_pooled_curves.tsv")
  )
  
  if (!file.exists(file)) {
    warning("Missing: ", file)
    return(
      ggplot() +
        theme_void() +
        labs(title = paste("B2 NR", trait, "\nMissing file"))
    )
  }
  
  d <- fread(file)
  
  d[, population := factor(
    population,
    levels = c("Wild", "Landrace", "Improved")
  )]
  
  ggplot(d, aes(x = anc_effect, fill = population, color = population)) +
    geom_density(
      alpha = 0.35,
      adjust = 1.2,
      linewidth = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = cols_b2, drop = FALSE) +
    scale_color_manual(values = cols_b2, drop = FALSE) +
    scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.2)) +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("B2 NR", trait),
      x = "Ancestral allele effect size",
      y = "Density",
      fill = "Population",
      color = "Population"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

#-----------------------------
# helper function: PosSel plot
#-----------------------------
make_pos_plot <- function(trait) {
  
  file <- file.path(
    BASE_DIR_POS,
    paste0("NR_", trait, "_PosSel_top1_pooled_curves.tsv")
  )
  
  if (!file.exists(file)) {
    warning("Missing: ", file)
    return(
      ggplot() +
        theme_void() +
        labs(title = paste("PosSel NR", trait, "\nMissing file"))
    )
  }
  
  d <- fread(file)
  
  d[, population := factor(
    population,
    levels = c("Landrace", "Improved")
  )]
  
  ggplot(d, aes(x = anc_effect, fill = population, color = population)) +
    geom_density(
      alpha = 0.40,
      adjust = 1.2,
      linewidth = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = cols_pos, drop = FALSE) +
    scale_color_manual(values = cols_pos, drop = FALSE) +
    scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.2)) +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("PosSel NR", trait),
      x = "Ancestral allele effect size",
      y = "Density",
      fill = "Population",
      color = "Population"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

#-----------------------------
# build plots
#-----------------------------
b2_plots  <- lapply(traits, make_b2_plot)
pos_plots <- lapply(traits, make_pos_plot)

#-----------------------------
# combine: 2 rows x 3 columns
#-----------------------------
combined_plot <- (b2_plots[[1]] + b2_plots[[2]] + b2_plots[[3]]) /
  (pos_plots[[1]] + pos_plots[[2]] + pos_plots[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#-----------------------------
# save one PDF
#-----------------------------
ggsave(
  filename = file.path(OUTDIR, "NR_B2_PosSel_density_2row_3col_zoomed.pdf"),
  plot = combined_plot,
  width = 14,
  height = 5,
  dpi = 600
)

# print
combined_plot