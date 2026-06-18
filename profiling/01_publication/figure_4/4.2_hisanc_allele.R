#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"

base_b2 <- file.path(repo_root, "data/4.GWAS_selection/A.balancing_selection")
base_ps <- file.path(repo_root, "data/4.GWAS_selection/B.positive_selection")
out_dir <- file.path(repo_root, "graphs/01_publication/4.GWAS_selection_integration")
out_pdf <- file.path(out_dir, "NR_B2_PosSel_histogram_2row_3col.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("Architecture", "Developmental", "Panicle")
trait_labels <- c(
  Architecture = "Architectural",
  Developmental = "Developmental",
  Panicle = "Panicle"
)

cols_b2_line <- c(
  Wild = "#1b9e77",
  Landrace = "#d95f02",
  Improved = "#7570b3"
)

cols_ps_line <- c(
  Domestication = "#E69F00",
  Improvement = "#009E73"
)

common_breaks <- seq(-0.3, 0.3, by = 0.02)

read_b2 <- function(trait) {
  file <- file.path(base_b2, paste0("NR_", trait, "_B2_top1_pooled_curves.tsv"))
  d <- fread(file)
  d[, trait := factor(
    trait_labels[[trait]],
    levels = c("Architectural", "Developmental", "Panicle")
  )]
  d[, population := factor(population, levels = c("Wild", "Landrace", "Improved"))]
  d
}

read_ps <- function(trait) {
  file <- file.path(base_ps, paste0("NR_", trait, "_PosSel_top1_raw_points.tsv"))
  d <- fread(file)
  d[, population := fifelse(population == "Landrace", "Domestication", "Improvement")]
  d[, trait := factor(
    trait_labels[[trait]],
    levels = c("Architectural", "Developmental", "Panicle")
  )]
  d[, population := factor(population, levels = c("Domestication", "Improvement"))]
  d
}

b2 <- rbindlist(lapply(traits, read_b2), fill = TRUE)
ps <- rbindlist(lapply(traits, read_ps), fill = TRUE)

make_row_plot <- function(d, fill_vals, line_vals, legend_inside) {
  ggplot(d, aes(x = anc_effect, fill = population, color = population)) +
    geom_histogram(
      aes(y = after_stat(density)),
      breaks = common_breaks,
      position = "identity",
      alpha = 0.35,
      linewidth = 0.35
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~trait, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    scale_color_manual(values = line_vals, drop = FALSE) +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    scale_x_continuous(breaks = c(-0.3, -0.1, 0.1, 0.3)) +
    labs(x = NULL, y = "Density", fill = NULL, color = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 14),
      legend.position = "inside",
      legend.position.inside = legend_inside,
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.margin = margin(5.5, 20, 5.5, 5.5)
    )
}

p_top <- make_row_plot(b2, cols_b2_line, cols_b2_line, c(0.08, 0.78))
p_bottom <- make_row_plot(ps, cols_ps_line, cols_ps_line, c(0.11, 0.78)) +
  labs(x = "Ancestral allele effect size")

combined_plot <- p_top / p_bottom + plot_layout(heights = c(1, 1))

ggsave(out_pdf, combined_plot, width = 14, height = 5, device = cairo_pdf)

message("Saved: ", out_pdf)
