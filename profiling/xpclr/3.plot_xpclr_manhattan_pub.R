library(data.table)
library(ggplot2)

repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"

chr_file <- file.path(repo_root, "data/2.2_positive_selec/sorg.chr_length_V5.txt")
dom_file <- file.path(
  repo_root,
  "data/xpclr/domestication_landrace_vs_wild/domestication_landrace_vs_wild.all_chr.tsv"
)
breed_file <- file.path(
  repo_root,
  "data/xpclr/breeding_improved_vs_landrace/breeding_improved_vs_landrace.all_chr.tsv"
)
out_dir <- file.path(repo_root, "graphs/01_publication/xpclr")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ch <- fread(chr_file)

prep <- function(path, label) {
  d <- fread(path)
  stopifnot(all(c("chrom", "start", "stop", "xpclr") %in% names(d)))

  d[, chrom := as.integer(chrom)]
  d[, pos := floor((start + stop) / 2)]
  d[, offset := ch$cumsum2[match(chrom, ch$Chr)]]
  d[, pos_cum := pos + offset]
  d[, pos_mb := pos_cum / 1e6]
  d[, chr_fill := factor(ifelse(chrom %% 2 == 1, "odd", "even"), levels = c("odd", "even"))]
  d[, comparison := label]
  d[is.finite(xpclr) & is.finite(pos_mb)]
}

dom <- prep(dom_file, "Landrace vs Wild")
breed <- prep(breed_file, "Improved vs Landrace")
plot_dt <- rbindlist(list(dom, breed), use.names = TRUE)

thresholds <- plot_dt[, .(
  top1pct_threshold = as.numeric(quantile(xpclr, 0.99, na.rm = TRUE))
), by = comparison]

x_ticks <- ch[, .(
  chrom = Chr,
  center_mb = (cumsum2 + Pos / 2) / 1e6
)]

p <- ggplot(plot_dt, aes(x = pos_mb, y = xpclr, color = chr_fill)) +
  geom_point(size = 0.28, alpha = 0.55, stroke = 0) +
  geom_hline(
    data = thresholds,
    aes(yintercept = top1pct_threshold),
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_grid(comparison ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = c(odd = "#00000099", even = "#C9C9C9CC")) +
  scale_x_continuous(
    breaks = x_ticks$center_mb,
    labels = x_ticks$chrom,
    limits = c(0, max(ch$cumsum1) / 1e6),
    expand = expansion(mult = c(0.002, 0.01))
  ) +
  labs(x = "Chromosome", y = "XP-CLR") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    strip.placement = "outside",
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(out_dir, "xpclr_manhattan_pub.pdf"),
  p,
  width = 14,
  height = 7,
  units = "in"
)

ggsave(
  file.path(out_dir, "xpclr_manhattan_pub.tiff"),
  p,
  width = 14,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

fwrite(
  thresholds,
  file.path(out_dir, "xpclr_manhattan_pub_top1pct_thresholds.txt"),
  sep = "\t"
)
