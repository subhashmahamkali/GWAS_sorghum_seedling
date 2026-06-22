#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

input_dir <- Sys.getenv(
  "ARCHIVE_5PCT_INPUT_DIR",
  unset = "/Users/subhashmahamkali/Downloads/Archive (1)"
)
out_dir <- Sys.getenv(
  "ARCHIVE_5PCT_OUT_DIR",
  unset = "/Users/subhashmahamkali/Documents/New project/archive_5pct_sweeps_delterious_plot"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

observed_file <- file.path(input_dir, "Sweeps_Deleterious_SNP_num_kb_observed.txt")
comparison_files <- c(
  "Landrace_vs_Wild_merged_sampling_delterious.txt",
  "Improved_vs_Landrace_merged_sampling_delterious.txt"
)
group_files <- c(
  "B2_sweeps_wild_sampling_delterious.txt",
  "B2_sweeps_landrace_sampling_delterious.txt",
  "B2_sweeps_modern_sampling_delterious.txt"
)

required_files <- c(
  observed_file,
  file.path(input_dir, comparison_files),
  file.path(input_dir, group_files)
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing input files:\n", paste0("  - ", missing_files, collapse = "\n"))
}

type_map <- data.table(
  type = c(
    "Landrace_vs_Wild_merged",
    "Improved_vs_Landrace_merged",
    "B2_sweeps_wild",
    "B2_sweeps_landrace",
    "B2_sweeps_modern"
  ),
  panel = c("c", "c", "d", "d", "d"),
  group = c("lr_vs_wd", "md_vs_lr", "wd", "lr", "md"),
  label = c("LR vs. WD", "MD vs. LR", "WD", "LR", "MD"),
  color = c("#1E22FF", "#FF69B4", "#B45CFF", "#FFAC21", "#5A9E2F")
)
color_map <- setNames(type_map$color, type_map$group)
label_map <- setNames(type_map$label, type_map$group)

to_stars <- function(p) {
  if (is.na(p)) "ns"
  else if (p <= 0.01) "**"
  else if (p <= 0.05) "*"
  else "ns"
}

read_sampling_file <- function(path) {
  dt <- fread(path)
  base_type <- sub("_sampling_delterious.txt$", "", basename(path))
  dt[, type := base_type]
  dt[, sample_id := sub("^.*-Sample", "", Type)]
  dt[, sample_id := as.integer(sample_id)]
  dt[, .(
    type,
    sample_id,
    n_del = Deleterious_SNP_num,
    total_kb = `length(kb)`,
    per_kb = `Deleterious_SNP_num/kb`
  )]
}

obs_dt <- fread(observed_file)[, .(
  type = Type,
  obs_n_del = Deleterious_SNP_num,
  obs_total_kb = `length(kb)`,
  obs_per_kb = `Deleterious_SNP_num/kb`
)]

perm_dt <- rbindlist(
  lapply(file.path(input_dir, c(comparison_files, group_files)), read_sampling_file),
  use.names = TRUE
)

perm_dt <- merge(perm_dt, type_map, by = "type", all.x = TRUE, sort = FALSE)
obs_dt <- merge(obs_dt, type_map, by = "type", all.x = TRUE, sort = FALSE)

if (any(is.na(perm_dt$group)) || any(is.na(obs_dt$group))) {
  stop("Failed to map one or more input `Type` values to plotting groups.")
}

summary_dt <- perm_dt[, .(
  n_perm = .N,
  perm_median = median(per_kb, na.rm = TRUE),
  perm_mean = mean(per_kb, na.rm = TRUE),
  perm_min = min(per_kb, na.rm = TRUE),
  perm_max = max(per_kb, na.rm = TRUE)
), by = .(type, panel, group, label, color)]

summary_dt <- merge(summary_dt, obs_dt, by = c("type", "panel", "group", "label", "color"), all.x = TRUE)
summary_dt <- merge(
  summary_dt,
  perm_dt[, .(n_perm_total = .N), by = type],
  by = "type",
  all.x = TRUE
)
summary_dt[, p_depletion := mapply(
  FUN = function(tt, obs_val, n_perm_total) {
    ((perm_dt[type == tt, sum(per_kb <= obs_val)]) + 1) / (n_perm_total + 1)
  },
  tt = type,
  obs_val = obs_per_kb,
  n_perm_total = n_perm_total
)]
summary_dt[, star_label := vapply(p_depletion, to_stars, character(1))]
summary_dt[, cutoff := factor("5%", levels = "5%")]

perm_dt[, cutoff := factor("5%", levels = "5%")]

panel_c_levels <- type_map[panel == "c", group]
panel_d_levels <- type_map[panel == "d", group]

build_panel <- function(panel_id, group_levels, panel_tag) {
  pdt <- copy(perm_dt[panel == panel_id])
  sdt <- copy(summary_dt[panel == panel_id])
  
  pdt[, group := factor(group, levels = group_levels)]
  sdt[, group := factor(group, levels = group_levels)]
  
  offsets <- if (length(group_levels) == 2L) c(-0.20, 0.20) else c(-0.24, 0, 0.24)
  violin_width <- if (length(group_levels) == 2L) 0.18 else 0.18
  box_width <- if (length(group_levels) == 2L) 0.03 else 0.035
  line_halfwidth <- if (length(group_levels) == 2L) 0.07 else 0.08
  offset_map <- data.table(
    group = factor(group_levels, levels = group_levels),
    x_plot = 1 + offsets
  )
  pdt <- merge(pdt, offset_map, by = "group", sort = FALSE)
  sdt <- merge(sdt, offset_map, by = "group", sort = FALSE)
  sdt[, x_start := x_plot - line_halfwidth]
  sdt[, x_end := x_plot + line_halfwidth]
  
  y_min <- min(c(pdt$per_kb, sdt$obs_per_kb), na.rm = TRUE)
  y_max <- max(c(pdt$per_kb, sdt$obs_per_kb), na.rm = TRUE)
  y_pad <- max((y_max - y_min) * 0.12, 0.06)
  star_pad <- max((y_max - y_min) * 0.06, 0.04)
  
  ggplot(pdt, aes(x = x_plot, y = per_kb, fill = group, color = group)) +
    geom_violin(width = violin_width, alpha = 1, linewidth = 0.75, trim = FALSE) +
    geom_boxplot(
      width = box_width,
      fill = "white",
      outlier.shape = NA,
      linewidth = 0.45
    ) +
    geom_segment(
      data = sdt,
      aes(x = x_start, xend = x_end, y = obs_per_kb, yend = obs_per_kb, color = group),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    geom_text(
      data = sdt,
      aes(x = x_plot, y = obs_per_kb + star_pad, label = star_label),
      inherit.aes = FALSE,
      color = "red",
      fontface = "bold",
      size = 6
    ) +
    scale_fill_manual(values = color_map[group_levels], labels = label_map[group_levels]) +
    scale_color_manual(values = color_map[group_levels], labels = label_map[group_levels]) +
    scale_x_continuous(
      breaks = 1,
      labels = "5%",
      limits = c(min(offset_map$x_plot) - 0.30, max(offset_map$x_plot) + 0.30)
    ) +
    coord_cartesian(ylim = c(max(0, y_min - y_pad * 0.35), y_max + y_pad)) +
    labs(x = NULL, y = "NO. of deleterious SNPs per kb", tag = panel_tag) +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.5),
      panel.grid.minor = element_line(color = "#F0F0F0", linewidth = 0.35),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.03, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_blank()
    )
}

panel_c <- build_panel("c", panel_c_levels, NULL)
panel_d <- build_panel("d", panel_d_levels, NULL)

combined_plot <- panel_c + panel_d + plot_layout(widths = c(1, 1.45))

summary_out <- summary_dt[, .(
  panel,
  label,
  type,
  n_perm,
  obs_per_kb,
  perm_median,
  perm_mean,
  perm_min,
  perm_max,
  p_depletion,
  star_label
)]

fwrite(
  summary_out,
  file.path(out_dir, "archive_5pct_sweeps_delterious_summary.tsv"),
  sep = "\t"
)

ggsave(
  filename = file.path(out_dir, "archive_5pct_sweeps_delterious_violin.pdf"),
  plot = combined_plot,
  width = 12.8,
  height = 5.7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "archive_5pct_sweeps_delterious_violin.png"),
  plot = combined_plot,
  width = 12.8,
  height = 5.7,
  units = "in",
  dpi = 350,
  bg = "white"
)

message("Saved plot and summary to: ", out_dir)
