#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------
# Inputs (from HCC pipeline)
# ------------------------------------------------------------------
in_zss <- "/mnt/nrdstor/jyanglab/subhash/datasets/3.results/negative/zss_neg_feature_exclusive_intergenic_for_thresholds.tsv"
in_bp  <- "/mnt/nrdstor/jyanglab/subhash/datasets/3.results/negative/genomic_fetures/feature_bp_exclusive_with_intergenic.tsv"
outdir <- "/mnt/nrdstor/jyanglab/subhash/datasets/3.results/negative/genomic_fetures/plots_1_5_10_with_intergenic"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

feature_levels <- c("exon", "intron", "upstream_2kb", "downstream_2kb", "intergenic")
cut_probs <- c("1%" = 0.01, "5%" = 0.05, "10%" = 0.10)

pal <- c(
  exon = "#4DAF4A",
  intron = "#377EB8",
  upstream_2kb = "#984EA3",
  downstream_2kb = "#E69F00",
  intergenic = "#8E8E8E"
)

z <- fread(in_zss, col.names = c("snpid", "zss", "feature"))
z <- z[is.finite(zss) & zss < 0 & feature %in% feature_levels]
z[, feature := factor(feature, levels = feature_levels)]

bp <- fread(in_bp, col.names = c("feature", "bp"))
bp <- bp[feature %in% feature_levels]
bp[, feature := factor(feature, levels = feature_levels)]

# thresholds on ALL negative SNPs (same style as earlier analyses)
cuts <- quantile(z$zss, probs = cut_probs, na.rm = TRUE, names = FALSE, type = 7)
names(cuts) <- names(cut_probs)

summ_list <- list()
plot_list <- list()

for (lab in names(cuts)) {
  cval <- cuts[[lab]]
  zz <- z[zss <= cval]

  s <- zz[, .(
    n_neg = .N,
    median_zss = median(zss),
    mean_zss = mean(zss),
    sum_neg_zss = sum(zss)
  ), by = feature][bp, on = "feature"]

  s[is.na(n_neg), `:=`(n_neg = 0, median_zss = NA_real_, mean_zss = NA_real_, sum_neg_zss = 0)]
  s[, cutoff := lab]
  s[, threshold_zss := cval]
  s[, neg_snp_per_mb := n_neg / (bp / 1e6)]
  s[, neg_snp_per_kb := n_neg / (bp / 1e3)]
  summ_list[[lab]] <- s

  # sample for plotting to keep file manageable
  ps <- zz[, .SD[sample(.N, min(.N, 120000))], by = feature]
  ps[, cutoff := lab]
  ps[, threshold_zss := cval]
  plot_list[[lab]] <- ps
}

summ <- rbindlist(summ_list, use.names = TRUE, fill = TRUE)
plot_dt <- rbindlist(plot_list, use.names = TRUE, fill = TRUE)

summ[, cutoff := factor(cutoff, levels = c("1%", "5%", "10%"))]
plot_dt[, cutoff := factor(cutoff, levels = c("1%", "5%", "10%"))]

fwrite(summ, file.path(outdir, "zss_feature_threshold_summary_with_intergenic.tsv"), sep = "\t")
fwrite(data.table(cutoff = names(cuts), zss_threshold = as.numeric(cuts)),
       file.path(outdir, "zss_thresholds_1_5_10.tsv"), sep = "\t")

# Plot 1: Advisor-requested normalized count
p1 <- ggplot(summ, aes(x = feature, y = neg_snp_per_mb, fill = feature)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.22) +
  facet_wrap(~ cutoff, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Negative SNP Density by Genomic Feature (normalized per Mb)",
    x = NULL,
    y = "Negative SNPs per Mb"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", angle = 18, hjust = 1)
  )

# Plot 2: ZSS distribution by feature at each threshold
p2 <- ggplot(plot_dt, aes(x = feature, y = zss, fill = feature)) +
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.2, alpha = 0.9) +
  geom_boxplot(width = 0.10, outlier.size = 0.12, fill = "white", color = "black", linewidth = 0.22) +
  facet_wrap(~ cutoff, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(
    title = "Distribution of Negative ZSS by Genomic Feature",
    x = NULL,
    y = "Negative ZSS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", angle = 18, hjust = 1)
  )

ggsave(file.path(outdir, "neg_snp_per_mb_by_feature_1_5_10_with_intergenic.pdf"),
       p1, width = 12.5, height = 4.6, units = "in", bg = "white", useDingbats = FALSE)
ggsave(file.path(outdir, "zss_distribution_by_feature_violin_1_5_10_with_intergenic.pdf"),
       p2, width = 12.5, height = 4.8, units = "in", bg = "white", useDingbats = FALSE)

message("Saved summary and plots to: ", outdir)
