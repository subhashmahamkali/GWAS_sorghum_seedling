#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(patchwork)
})

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"

three_species_dir <- file.path(
  repo_root,
  "graphs/01_publication/3.genetic_load/genome_level_1_5_10_violin_custom_palette"
)
sap_dir <- file.path(
  repo_root,
  "graphs/01_publication/3.genetic_load/genome_level_top1_four_species_violin"
)
out_dir <- file.path(
  repo_root,
  "graphs/01_publication/3.genetic_load/combined_A4_panels"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("wild", "landrace", "improved")

three_cols <- c(
  wild = "#4DAF4A",
  landrace = "#377EB8",
  improved = "#984EA3"
)
sap_cols <- c(
  "Durra" = "#F28E2B",
  "Kafir" = "#DCA73C",
  "Caudatum" = "#A56CC1",
  "Guinea" = "#77BCE8",
  "Milo/Durra-Bicolor" = "#E88AC6",
  "Bicolor" = "#B9B9B9"
)

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#7F7F7F40", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(3, 4, 3, 4)
  )

comparisons <- list(
  c("wild", "landrace"),
  c("wild", "improved"),
  c("landrace", "improved")
)

make_three_group_panel <- function(dt, species_label, cutoff_label) {
  d <- copy(dt[cutoff_lab == cutoff_label])
  d[, group := factor(as.character(group), levels = group_levels)]

  ggplot(d, aes(x = group, y = Abs_genetic_load, fill = group)) +
    geom_violin(width = 0.9, trim = FALSE, color = NA, linewidth = 0, alpha = 0.95) +
    geom_boxplot(width = 0.16, outlier.size = 0.35, fill = "white", color = "black", linewidth = 0.3) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 3.4,
      tip.length = 0.01
    ) +
    scale_fill_manual(values = three_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(title = paste(species_label, cutoff_label), y = "Absolute genetic load") +
    base_theme
}

make_sap_panel <- function(cutoff_label) {
  sap_file <- file.path(sap_dir, paste0("SAP_top", gsub("%", "", cutoff_label), "_by_sample.tsv"))
  let_file <- file.path(sap_dir, paste0("SAP_top", gsub("%", "", cutoff_label), "_compact_letters.tsv"))

  d <- fread(sap_file)
  d <- d[cutoff == cutoff_label & !is.na(Race) & !is.na(Abs_genetic_load)]
  ord <- d[, .(med = median(Abs_genetic_load, na.rm = TRUE)), by = Race][order(-med), as.character(Race)]
  d[, Race := factor(as.character(Race), levels = ord)]

  letters_dt <- fread(let_file)
  letters_dt[, Race := factor(as.character(Race), levels = ord)]
  y_rng <- range(d$Abs_genetic_load, na.rm = TRUE)
  y_pad <- ifelse(diff(y_rng) > 0, 0.05 * diff(y_rng), 0.05 * abs(y_rng[2]))
  y_pos <- d[, .(y = max(Abs_genetic_load, na.rm = TRUE) + y_pad), by = Race]
  letters_dt <- merge(letters_dt, y_pos, by = "Race", all.x = TRUE)

  ggplot(d, aes(x = Race, y = Abs_genetic_load, fill = Race)) +
    geom_violin(width = 0.9, trim = FALSE, color = NA, linewidth = 0, alpha = 0.95) +
    geom_boxplot(width = 0.16, outlier.size = 0.35, fill = "white", color = "black", linewidth = 0.3) +
    geom_text(
      data = letters_dt,
      aes(x = Race, y = y, label = cld),
      inherit.aes = FALSE,
      size = 3.4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = sap_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(title = paste("SAP", cutoff_label), y = "Absolute genetic load") +
    base_theme +
    theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1, face = "bold"))
}

make_sap_panel_custom_palette <- function(cutoff_label, custom_cols) {
  sap_file <- file.path(sap_dir, paste0("SAP_top", gsub("%", "", cutoff_label), "_by_sample.tsv"))
  let_file <- file.path(sap_dir, paste0("SAP_top", gsub("%", "", cutoff_label), "_compact_letters.tsv"))

  d <- fread(sap_file)
  d <- d[cutoff == cutoff_label & !is.na(Race) & !is.na(Abs_genetic_load)]
  ord <- d[, .(med = median(Abs_genetic_load, na.rm = TRUE)), by = Race][order(-med), as.character(Race)]
  d[, Race := factor(as.character(Race), levels = ord)]

  letters_dt <- fread(let_file)
  letters_dt[, Race := factor(as.character(Race), levels = ord)]
  y_rng <- range(d$Abs_genetic_load, na.rm = TRUE)
  y_pad <- ifelse(diff(y_rng) > 0, 0.05 * diff(y_rng), 0.05 * abs(y_rng[2]))
  y_pos <- d[, .(y = max(Abs_genetic_load, na.rm = TRUE) + y_pad), by = Race]
  letters_dt <- merge(letters_dt, y_pos, by = "Race", all.x = TRUE)

  ggplot(d, aes(x = Race, y = Abs_genetic_load, fill = Race)) +
    geom_violin(width = 0.9, trim = FALSE, color = NA, linewidth = 0, alpha = 0.95) +
    geom_boxplot(width = 0.16, outlier.size = 0.35, fill = "white", color = "black", linewidth = 0.3) +
    geom_text(
      data = letters_dt,
      aes(x = Race, y = y, label = cld),
      inherit.aes = FALSE,
      size = 3.4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = custom_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(title = paste("SAP", cutoff_label), y = "Absolute genetic load") +
    base_theme +
    theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1, face = "bold"))
}

sorghum <- fread(file.path(three_species_dir, "sorghum/sorghum_genome_level_by_sample_1_5_10.tsv"))
maize <- fread(file.path(three_species_dir, "maize/maize_genome_level_by_sample_1_5_10.tsv"))
sunflower <- fread(file.path(three_species_dir, "sunflower/sunflower_genome_level_by_sample_1_5_10.tsv"))

panel_list <- list(
  make_three_group_panel(sorghum, "Sorghum", "5%"),
  make_three_group_panel(maize, "Maize", "5%"),
  make_three_group_panel(sunflower, "Sunflower", "5%"),
  make_sap_panel("1%"),
  make_sap_panel("5%"),
  make_sap_panel("10%")
)

combined <- wrap_plots(panel_list, ncol = 3, byrow = TRUE) +
  plot_annotation(
    title = "Genome-Wide Absolute Genetic Load",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

pdf_out <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_A4.pdf")
png_out <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_A4.png")
pdf_out_60 <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_A4_60pct.pdf")
png_out_60 <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_A4_60pct.png")
pdf_out_21x20 <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_21x20cm.pdf")
png_out_21x20 <- file.path(out_dir, "genetic_load_species5_then_SAP_1_5_10_21x20cm.png")

w_in <- 8.27
h_in <- 11.69
scale_60 <- 0.60

ggsave(pdf_out, combined, width = w_in, height = h_in, units = "in", bg = "white", useDingbats = FALSE)
ggsave(png_out, combined, width = w_in, height = h_in, units = "in", dpi = 350, bg = "white")
ggsave(pdf_out_60, combined, width = w_in * scale_60, height = h_in * scale_60, units = "in", bg = "white", useDingbats = FALSE)
ggsave(png_out_60, combined, width = w_in * scale_60, height = h_in * scale_60, units = "in", dpi = 350, bg = "white")
ggsave(pdf_out_21x20, combined, width = 21, height = 20, units = "cm", bg = "white", useDingbats = FALSE)
ggsave(png_out_21x20, combined, width = 21, height = 20, units = "cm", dpi = 350, bg = "white")

# Panel B only (SAP 1/5/10) with requested race palette
requested_race_pal <- c(
  "durra" = "#FF7F00",
  "kafir" = "#E69F00",
  "caudatum" = "#984EA3",
  "milo/durra-bicolor" = "#F781BF",
  "mixed/bicolor" = "#A6A6A6",
  "guinea" = "#56B4E9"
)
requested_sap_cols <- c(
  "Durra" = requested_race_pal[["durra"]],
  "Kafir" = requested_race_pal[["kafir"]],
  "Caudatum" = requested_race_pal[["caudatum"]],
  "Milo/Durra-Bicolor" = requested_race_pal[["milo/durra-bicolor"]],
  "Bicolor" = requested_race_pal[["mixed/bicolor"]],
  "Guinea" = requested_race_pal[["guinea"]]
)

panel_b_only <- wrap_plots(
  list(
    make_sap_panel_custom_palette("1%", requested_sap_cols),
    make_sap_panel_custom_palette("5%", requested_sap_cols),
    make_sap_panel_custom_palette("10%", requested_sap_cols)
  ),
  ncol = 3,
  byrow = TRUE
) +
  plot_annotation(
    title = "SAP Absolute Genetic Load",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

panel_b_pdf <- file.path(out_dir, "genetic_load_panelB_SAP_1_5_10_race_palette_21cm.pdf")
panel_b_png <- file.path(out_dir, "genetic_load_panelB_SAP_1_5_10_race_palette_21cm.png")
ggsave(panel_b_pdf, panel_b_only, width = 21, height = 7.2, units = "cm", bg = "white", useDingbats = FALSE)
ggsave(panel_b_png, panel_b_only, width = 21, height = 7.2, units = "cm", dpi = 350, bg = "white")

message("Saved: ", pdf_out)
message("Saved: ", panel_b_pdf)
