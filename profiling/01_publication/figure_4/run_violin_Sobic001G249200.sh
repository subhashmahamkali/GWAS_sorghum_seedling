#!/bin/bash
# ============================================================================
# Step 1: Extract lead SNP genotype from PLINK
# Step 2: Run R to make violin plot for Sobic.001G249200
# ============================================================================
# Usage: bash run_violin_Sobic001G249200.sh
# ============================================================================

set -euo pipefail

# --- Paths ------------------------------------------------------------------
PLINK_PREFIX="/mnt/nrdstor/jyanglab/subhash/datasets/1.SAP_SNPs/2.filtered_vcf/1.het_0.1_SAP"
OUT_DIR="/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/zoom_inplot"
PHENO_FILE="/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/1.BLUEs_SAP_2020_2021.csv"

# Lead SNP info for Sobic.001G249200 / third leaf length NR
SNP_ID="Chr01_27091166"
OUT_PREFIX="${OUT_DIR}/leadSNP_Sobic001G249200"
awk '$1=1 && $2==27091166' /work/jyanglab/subhash/sorgsd/bal_s/ans_freq/Sorghum_ancestral_allele_V3.1.txt

mkdir -p "${OUT_DIR}"

# --- Step 1: Extract genotype with PLINK ------------------------------------
echo "=== Extracting genotype at ${SNP_ID} ==="

plink \
  --bfile "${PLINK_PREFIX}" \
  --snp "${SNP_ID}" \
  --recode A \
  --out "${OUT_PREFIX}"

echo "=== Genotype extracted: ${OUT_PREFIX}.raw ==="
echo ""

# --- Step 2: Run R violin plot ----------------------------------------------
echo "=== Running R violin plot ==="

Rscript - <<'RSCRIPT'

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# --- Paths ------------------------------------------------------------------
out_dir    <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/zoom_inplot"
geno_file  <- file.path(out_dir, "leadSNP_Sobic001G249200.raw")
pheno_file <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/1.BLUEs_SAP_2020_2021.csv"
bim_file   <- "/mnt/nrdstor/jyanglab/subhash/datasets/1.SAP_SNPs/2.filtered_vcf/1.het_0.1_SAP.bim"

snp_id      <- "Chr01_27091166"
gene_name   <- "Sobic.001G249200"
gene_desc   <- "KH domain"
trait_col   <- "X2020.thirdLeafLength.NR"
trait_label <- "Third leaf length (NR)"

# --- Read genotype ----------------------------------------------------------
cat("Reading genotype data...\n")
geno <- read.table(geno_file, header = TRUE, stringsAsFactors = FALSE)

dosage_col <- grep(paste0("^", snp_id), names(geno), value = TRUE)
if (length(dosage_col) == 0) {
  stop("Lead SNP dosage column not found in .raw file: ", snp_id)
}
dosage_col <- dosage_col[1]
cat(sprintf("  Dosage column: %s\n", dosage_col))

geno_df <- data.frame(
  sample = gsub("\\s+", "", geno$IID),
  dosage = geno[[dosage_col]],
  stringsAsFactors = FALSE
)

# --- Read phenotype ---------------------------------------------------------
cat("Reading phenotype data...\n")
pheno <- read.csv(pheno_file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE)
pheno$genotype <- gsub("\\s+", "", pheno$genotype)

if (!trait_col %in% names(pheno)) {
  trait_matches <- grep("thirdLeafLength.*NR$", names(pheno), value = TRUE)
  if (length(trait_matches) == 0) {
    stop("Phenotype column not found. Expected: ", trait_col)
  }
  trait_col <- tail(trait_matches, 1)
}
cat(sprintf("  Phenotype column: %s\n", trait_col))

pheno_df <- data.frame(
  sample = pheno$genotype,
  trait  = as.numeric(pheno[[trait_col]]),
  stringsAsFactors = FALSE
)

# --- Merge ------------------------------------------------------------------
merged <- merge(geno_df, pheno_df, by = "sample")
merged <- merged[complete.cases(merged), ]
cat(sprintf("  Merged samples: %d\n", nrow(merged)))

# --- Assign genotype labels -------------------------------------------------
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
snp_row <- bim[bim$V2 == snp_id, ]

if (nrow(snp_row) > 0) {
  # In .bim: V5 = A1 counted allele in --recode A, V6 = other allele
  a1 <- snp_row$V5
  a2 <- snp_row$V6
  cat(sprintf("  A1 counted allele: %s, A2 other allele: %s\n", a1, a2))

  merged$geno_label <- ifelse(
    merged$dosage == 0, paste0(a2, a2),
    ifelse(merged$dosage == 1, paste0(a1, a2),
           ifelse(merged$dosage == 2, paste0(a1, a1), NA))
  )
} else {
  warning("SNP not found in .bim; using generic genotype labels.")
  merged$geno_label <- ifelse(
    merged$dosage == 0, "REF/REF",
    ifelse(merged$dosage == 1, "REF/ALT",
           ifelse(merged$dosage == 2, "ALT/ALT", NA))
  )
}

# --- Drop heterozygotes for 2-class comparison ------------------------------
merged <- merged[merged$dosage != 1, ]
merged <- merged[complete.cases(merged), ]
cat(sprintf("  Homozygous samples only: %d\n", nrow(merged)))

ref_label <- unique(merged$geno_label[merged$dosage == 0])
alt_label <- unique(merged$geno_label[merged$dosage == 2])

if (length(ref_label) == 0 || length(alt_label) == 0) {
  stop("Need both homozygous genotype classes after filtering heterozygotes.")
}

ref_label <- ref_label[1]
alt_label <- alt_label[1]
merged$geno_label <- factor(merged$geno_label, levels = c(ref_label, alt_label))

# --- Statistics -------------------------------------------------------------
grp0 <- merged$trait[merged$dosage == 0]
grp2 <- merged$trait[merged$dosage == 2]

if (length(grp0) < 2 || length(grp2) < 2) {
  warning("One homozygous group has fewer than 2 samples; Wilcoxon test may be weak.")
}

wt <- wilcox.test(grp0, grp2)
pval <- wt$p.value
sig_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"

cat(sprintf("  Wilcoxon p-value: %.4e (%s)\n", pval, sig_label))
cat(sprintf("  n(%s) = %d, n(%s) = %d\n", ref_label, length(grp0), alt_label, length(grp2)))
cat(sprintf("  mean(%s) = %.3f, mean(%s) = %.3f\n", ref_label, mean(grp0), alt_label, mean(grp2)))

n_tab <- table(merged$geno_label)
n_labels <- paste0("n=", as.integer(n_tab[levels(merged$geno_label)]))

# --- Plot settings ----------------------------------------------------------
fill_cols <- c("#66C2A5", "#8DA0CB")
names(fill_cols) <- c(ref_label, alt_label)

ymax <- max(merged$trait, na.rm = TRUE)
ymin <- min(merged$trait, na.rm = TRUE)
yrange <- ymax - ymin
if (yrange == 0) yrange <- 1

text_y    <- ymax + yrange * 0.05
bracket_y <- ymax + yrange * 0.13
sig_y     <- ymax + yrange * 0.17
plot_top  <- ymax + yrange * 0.23

# --- Build publication-style violin plot ------------------------------------
cat("Generating violin plot...\n")

p <- ggplot(merged, aes(x = geno_label, y = trait, fill = geno_label)) +
  geom_violin(
    trim = TRUE,
    width = 0.82,
    color = "grey25",
    linewidth = 0.45,
    alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = fill_cols) +
  annotate(
    "text",
    x = 1:2,
    y = rep(text_y, 2),
    label = n_labels,
    size = 3.5
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = bracket_y, yend = bracket_y,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = bracket_y, yend = bracket_y - yrange * 0.02,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = bracket_y, yend = bracket_y - yrange * 0.02,
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = sig_y,
    label = sig_label,
    size = 4.8,
    fontface = "bold"
  ) +
  coord_cartesian(ylim = c(ymin, plot_top), clip = "off") +
  labs(
    x = "Genotype",
    y = trait_label,
    title = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 13, face = "bold", colour = "black"),
    axis.text.x  = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 11, colour = "black"),
    axis.line    = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    panel.grid   = element_blank(),
    plot.margin  = margin(10, 10, 10, 10)
  )

# --- Save outputs -----------------------------------------------------------
out_pdf  <- file.path(out_dir, "Sobic001G249200_violin_genotype_effect_MBEstyle_nojitter.pdf")
out_tiff <- file.path(out_dir, "Sobic001G249200_violin_genotype_effect_MBEstyle_nojitter.tiff")

ggsave(out_pdf, p, width = 4.0, height = 4.8, device = cairo_pdf)
ggsave(out_tiff, p, width = 4.0, height = 4.8, dpi = 600, compression = "lzw")

cat(sprintf("\nDone! Files saved to:\n  %s\n  %s\n", out_pdf, out_tiff))

RSCRIPT

echo "=== All done! ==="
