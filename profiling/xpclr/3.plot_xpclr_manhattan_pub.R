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
setnames(ch, c("Chr", "Pos", "cumsum1", "cumsum2"))

chr_centers <- ch[, .(
  chrom = Chr,
  center_mb = (cumsum2 + Pos / 2) / 1e6
)]

chr_boundaries <- data.table(
  x = c(ch$cumsum2[-1], tail(ch$cumsum1, 1)) / 1e6
)

prep <- function(path, label) {
  d <- fread(path)
  stopifnot(all(c("chrom", "start", "stop", "xpclr") %in% names(d)))

  d[, chrom := as.integer(chrom)]
  d[, midpoint := floor((start + stop) / 2)]
  d[, offset := ch$cumsum2[match(chrom, ch$Chr)]]
  d[, genome_bp := midpoint + offset]
  d[, genome_mb := genome_bp / 1e6]
  d[, chr_fill := factor(ifelse(chrom %% 2 == 1, "odd", "even"), levels = c("odd", "even"))]
  d[, comparison := label]

  d[is.finite(xpclr) & is.finite(genome_mb)]
}

dom <- prep(dom_file, "Landrace vs Wild")
breed <- prep(breed_file, "Improved vs Landrace")

plot_dt <- rbindlist(list(dom, breed), use.names = TRUE)
plot_dt[, comparison := factor(comparison, levels = c("Landrace vs Wild", "Improved vs Landrace"))]

thresholds <- plot_dt[, .(
  top1pct_threshold = as.numeric(quantile(xpclr, 0.99, na.rm = TRUE))
), by = comparison]

max_x <- max(ch$cumsum1) / 1e6

p <- ggplot(plot_dt, aes(x = genome_mb, y = xpclr, color = chr_fill)) +
  geom_vline(
    data = chr_boundaries,
    aes(xintercept = x),
    inherit.aes = FALSE,
    color = "#DDDDDD",
    linewidth = 0.35
  ) +
  geom_point(size = 0.22, alpha = 0.58, stroke = 0) +
  geom_hline(
    data = thresholds,
    aes(yintercept = top1pct_threshold),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 0.65
  ) +
  facet_grid(comparison ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = c(odd = "#00000080", even = "#C7C7C7CC")) +
  scale_x_continuous(
    breaks = chr_centers$center_mb,
    labels = chr_centers$chrom,
    limits = c(0, max_x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Chromosome", y = "XP-CLR") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.8, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 13, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 6)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.line.x = element_line(linewidth = 0.6, color = "black"),
    axis.line.y = element_line(linewidth = 0.6, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.5, color = "black"),
    panel.border = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  file.path(out_dir, "xpclr_manhattan_pub.pdf"),
  p,
  width = 14,
  height = 7,
  units = "in",
  bg = "white"
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
