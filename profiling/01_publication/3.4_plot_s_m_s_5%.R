#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(patchwork)
})

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"

out_dir <- file.path(
  repo_root,
  "graphs/01_publication/3.genetic_load/sorghum_maize_sunflower_5pct"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("wild", "landrace", "improved")
group_cols <- c(
  wild = "#4DAF4A",
  landrace = "#377EB8",
  improved = "#984EA3"
)

comparisons <- list(
  c("wild", "landrace"),
  c("wild", "improved"),
  c("landrace", "improved")
)

to_cutoff_num <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  u <- sort(unique(na.omit(y)))
  if (length(u) && all(u %in% c(1, 5, 10, 25))) y <- y / 100
  y
}

prepare_species_5pct <- function(files, group_col = "sub", modern_to_improved = FALSE) {
  dt <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  dt[, cutoff_num := to_cutoff_num(cutoff)]
  dt <- dt[
    cutoff_num == 0.05 &
      !is.na(sample) &
      !is.na(get(group_col)) &
      !is.na(Genetic_load)
  ]

  dt[, group := as.character(get(group_col))]
  if (modern_to_improved) dt[group == "modern", group := "improved"]
  dt <- dt[group %in% group_levels]

  out <- dt[, .(
    Genetic_load_sum = sum(Genetic_load, na.rm = TRUE),
    N_del_sum = sum(N_del, na.rm = TRUE),
    n_chr_rows = .N
  ), by = .(sample, group)]

  out[, Abs_genetic_load := abs(Genetic_load_sum)]
  out[, group := factor(group, levels = group_levels)]
  out[]
}

pairwise_wilcox <- function(d) {
  lev <- as.character(group_levels)
  cmb <- t(combn(lev, 2))
  out <- vector("list", nrow(cmb))
  for (i in seq_len(nrow(cmb))) {
    g1 <- cmb[i, 1]
    g2 <- cmb[i, 2]
    x <- d[group == g1, Abs_genetic_load]
    y <- d[group == g2, Abs_genetic_load]
    wt <- suppressWarnings(wilcox.test(x, y))
    out[[i]] <- data.table(group1 = g1, group2 = g2, p_value = wt$p.value)
  }
  pw <- rbindlist(out)
  pw[, p_adj_bh := p.adjust(p_value, method = "BH")]
  pw[order(p_adj_bh)]
}

make_panel <- function(d, title_txt) {
  ggplot(d, aes(x = group, y = Abs_genetic_load, fill = group)) +
    geom_violin(width = 0.9, trim = FALSE, color = "black", linewidth = 0.35, alpha = 0.95) +
    geom_boxplot(width = 0.16, outlier.size = 0.35, fill = "white", color = "black", linewidth = 0.3) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      tip.length = 0.01,
      size = 4
    ) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(title = title_txt, x = NULL, y = "Absolute genetic load") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

sorghum_files <- Sys.glob(file.path(repo_root, "data/3.1_genetic_load_sorghum/raw_*_miss50_chr*.csv"))
maize_files <- Sys.glob(file.path(repo_root, "data/3.3_genetic_load_maize/raw_*_chr*.csv"))
sun_files <- Sys.glob(file.path(repo_root, "data/3.2_genetic_load_sunflower/raw_*_miss50_chr*.csv"))

if (!length(sorghum_files)) stop("No sorghum files found.")
if (!length(maize_files)) stop("No maize files found.")
if (!length(sun_files)) stop("No sunflower files found.")

sorghum_5 <- prepare_species_5pct(sorghum_files, group_col = "sub", modern_to_improved = FALSE)
maize_5 <- prepare_species_5pct(maize_files, group_col = "sub", modern_to_improved = FALSE)
sun_5 <- prepare_species_5pct(sun_files, group_col = "sub", modern_to_improved = TRUE)

fwrite(sorghum_5, file.path(out_dir, "sorghum_5pct_genome_level_by_sample.tsv"), sep = "\t")
fwrite(maize_5, file.path(out_dir, "maize_5pct_genome_level_by_sample.tsv"), sep = "\t")
fwrite(sun_5, file.path(out_dir, "sunflower_5pct_genome_level_by_sample.tsv"), sep = "\t")

fwrite(pairwise_wilcox(sorghum_5), file.path(out_dir, "sorghum_5pct_pairwise_wilcox.tsv"), sep = "\t")
fwrite(pairwise_wilcox(maize_5), file.path(out_dir, "maize_5pct_pairwise_wilcox.tsv"), sep = "\t")
fwrite(pairwise_wilcox(sun_5), file.path(out_dir, "sunflower_5pct_pairwise_wilcox.tsv"), sep = "\t")

summary_dt <- rbindlist(list(
  sorghum_5[, .(
    species = "sorghum", group, n = .N,
    median_abs = median(Abs_genetic_load, na.rm = TRUE),
    mean_abs = mean(Abs_genetic_load, na.rm = TRUE),
    min_abs = min(Abs_genetic_load, na.rm = TRUE),
    max_abs = max(Abs_genetic_load, na.rm = TRUE)
  ), by = group],
  maize_5[, .(
    species = "maize", group, n = .N,
    median_abs = median(Abs_genetic_load, na.rm = TRUE),
    mean_abs = mean(Abs_genetic_load, na.rm = TRUE),
    min_abs = min(Abs_genetic_load, na.rm = TRUE),
    max_abs = max(Abs_genetic_load, na.rm = TRUE)
  ), by = group],
  sun_5[, .(
    species = "sunflower", group, n = .N,
    median_abs = median(Abs_genetic_load, na.rm = TRUE),
    mean_abs = mean(Abs_genetic_load, na.rm = TRUE),
    min_abs = min(Abs_genetic_load, na.rm = TRUE),
    max_abs = max(Abs_genetic_load, na.rm = TRUE)
  ), by = group]
), use.names = TRUE)
fwrite(summary_dt, file.path(out_dir, "s_m_s_5pct_summary.tsv"), sep = "\t")

p_sorghum <- make_panel(sorghum_5, "Sorghum 5%")
p_maize <- make_panel(maize_5, "Maize 5%")
p_sun <- make_panel(sun_5, "Sunflower 5%")

combined <- (p_sorghum | p_maize | p_sun) +
  plot_annotation(
    title = "Genome-Wide Absolute Genetic Load (Top 5%)",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(
  file.path(out_dir, "s_m_s_5pct_3panel_violin.pdf"),
  combined, width = 13.5, height = 5.1, units = "in",
  bg = "white", useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "s_m_s_5pct_3panel_violin.png"),
  combined, width = 13.5, height = 5.1, units = "in",
  dpi = 350, bg = "white"
)

message("Saved 5% sorghum/maize/sunflower outputs to: ", out_dir)
