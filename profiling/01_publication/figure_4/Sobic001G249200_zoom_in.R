#!/usr/bin/env Rscript
# ============================================================================
# Zoom-in plot: GWAS (2020 thirdLeafLength NR) + FST (landrace vs improved)
# around Sobic.001G249200
# ============================================================================
# Output: cairo_pdf to largedata/zoom_inplot/
# ============================================================================

# --- Paths (edit if needed) -------------------------------------------------
gwas_file <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/All_GWAS_sorg.txt"
fst_file  <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/figure_2/landrace_vs_improved.windowed.weir.fst"
out_dir   <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/zoom_inplot/"

# --- Gene & GWAS info -------------------------------------------------------
gene_name  <- "Sobic.001G249200"
gene_desc  <- "KH domain"
gene_start <- 27057733
gene_end   <- 27067348
gene_mid   <- (gene_start + gene_end) / 2

chr_num    <- 1
lead_snp   <- 27091166
lead_p     <- 8.373712e-07
gwas_trait <- "X2020.thirdLeafLength.NR"

# Zoom window: 500 kb on each side of gene.
window_pad   <- 500000
region_start <- gene_start - window_pad
region_end   <- gene_end + window_pad

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Read & subset GWAS data
# ============================================================================
cat("Reading GWAS data...\n")
gwas <- read.delim(gwas_file, header = TRUE, stringsAsFactors = FALSE)

gwas_sub <- gwas[gwas$Trait == gwas_trait &
                   gwas$CHROM == chr_num &
                   gwas$POS >= region_start &
                   gwas$POS <= region_end, ]

if (nrow(gwas_sub) == 0) {
  stop("No GWAS rows found. Check gwas_trait, chr_num, and region.")
}

gwas_sub$logP <- -log10(gwas_sub$P.value)
gwas_sub$pos_mb <- gwas_sub$POS / 1e6
cat(sprintf("  GWAS SNPs in window: %d\n", nrow(gwas_sub)))

# ============================================================================
# 2. Read & subset FST data
# ============================================================================
cat("Reading FST data...\n")
fst <- read.delim(fst_file, header = TRUE, stringsAsFactors = FALSE)

fst_sub <- fst[fst$CHROM == chr_num &
                 fst$BIN_START >= region_start &
                 fst$BIN_END <= region_end, ]

if (nrow(fst_sub) == 0) {
  stop("No FST windows found. Check fst_file, chr_num, and region.")
}

fst_sub$pos_mid <- (fst_sub$BIN_START + fst_sub$BIN_END) / 2
fst_sub$pos_mb <- fst_sub$pos_mid / 1e6
fst_sub$WEIGHTED_FST[fst_sub$WEIGHTED_FST < 0] <- 0
cat(sprintf("  FST windows in region: %d\n", nrow(fst_sub)))

# ============================================================================
# 3. Thresholds
# ============================================================================
gwas_thresh <- 5

cat("Computing genome-wide FST threshold (top 1%)...\n")
fst_all <- fst$WEIGHTED_FST
fst_all[fst_all < 0] <- 0
fst_thresh <- quantile(fst_all, 0.99, na.rm = TRUE)
cat(sprintf("  FST 99th percentile threshold: %.4f\n", fst_thresh))

# ============================================================================
# 4. Plot
# ============================================================================
cat("Generating plot...\n")
out_file <- file.path(out_dir, "Sobic001G249200_chr1_zoom_in_GWAS_FST.pdf")

cairo_pdf(out_file, width = 7, height = 7)

par(mfrow = c(2, 1),
    mar = c(1.5, 5, 2.5, 1),
    oma = c(3, 0, 2, 0),
    mgp = c(3, 0.7, 0))

xlim <- c(region_start / 1e6, region_end / 1e6)

# ========== Panel A: GWAS ===================================================
plot(gwas_sub$pos_mb, gwas_sub$logP,
     pch = 19, cex = 0.8,
     col = adjustcolor("#6A5ACD", alpha.f = 0.7),
     xlim = xlim,
     ylim = c(0, max(c(gwas_sub$logP, gwas_thresh + 1), na.rm = TRUE)),
     xlab = "",
     ylab = expression(-log[10](italic(P))),
     xaxt = "n",
     las = 1,
     cex.lab = 1.3, cex.axis = 1.1)

rect(gene_start / 1e6, par("usr")[3],
     gene_end / 1e6, par("usr")[4],
     col = adjustcolor("#DAA520", alpha.f = 0.25), border = NA)

points(gwas_sub$pos_mb, gwas_sub$logP,
       pch = 19, cex = 0.8,
       col = adjustcolor("#6A5ACD", alpha.f = 0.7))

abline(h = gwas_thresh, col = "red", lty = 2, lwd = 1.5)

lead_snp_mb <- lead_snp / 1e6
lead_snp_row <- gwas_sub[which.min(abs(gwas_sub$POS - lead_snp)), , drop = FALSE]
if (nrow(lead_snp_row) > 0) {
  points(lead_snp_row$pos_mb, lead_snp_row$logP,
         pch = 18, cex = 2.2, col = "#33A02C")
  arrows(x0 = lead_snp_mb + 0.03, y0 = lead_snp_row$logP + 0.6,
         x1 = lead_snp_mb, y1 = lead_snp_row$logP + 0.15,
         length = 0.08, lwd = 1.5, col = "#33A02C")
  text(lead_snp_mb + 0.035, lead_snp_row$logP + 0.7,
       labels = paste0("Lead SNP\np = ", format(lead_p, scientific = TRUE)),
       cex = 0.6, col = "#33A02C", adj = 0, font = 2)
}

text(gene_mid / 1e6, par("usr")[4] * 0.92,
     labels = gene_name, cex = 0.75, font = 4, col = "#33A02C")

mtext("A", side = 3, line = 0.5, at = par("usr")[1], adj = 0, font = 2, cex = 1.4)
axis(1, labels = FALSE)

# ========== Panel B: FST ====================================================
par(mar = c(4, 5, 1, 1))

plot(fst_sub$pos_mb, fst_sub$WEIGHTED_FST,
     type = "l", lwd = 1.2,
     col = "black",
     xlim = xlim,
     ylim = c(0, max(c(fst_sub$WEIGHTED_FST, fst_thresh * 1.2), na.rm = TRUE)),
     xlab = "",
     ylab = expression(F[ST]),
     las = 1,
     cex.lab = 1.3, cex.axis = 1.1)

rect(gene_start / 1e6, par("usr")[3],
     gene_end / 1e6, par("usr")[4],
     col = adjustcolor("#DAA520", alpha.f = 0.25), border = NA)

polygon(c(fst_sub$pos_mb[1], fst_sub$pos_mb, tail(fst_sub$pos_mb, 1)),
        c(0, fst_sub$WEIGHTED_FST, 0),
        col = adjustcolor("grey30", alpha.f = 0.15), border = NA)
lines(fst_sub$pos_mb, fst_sub$WEIGHTED_FST, lwd = 1.2, col = "black")

abline(h = fst_thresh, col = "red", lty = 2, lwd = 1.5)

mtext("B", side = 3, line = 0.3, at = par("usr")[1], adj = 0, font = 2, cex = 1.4)
mtext("Chr1 (Mb)", side = 1, line = 2.5, cex = 1.2)

mtext(paste0("Third leaf length (NR) - ", gene_desc, " (", gene_name, ")"),
      side = 3, outer = TRUE, line = 0.3, cex = 1.1, font = 2)

dev.off()

cat(sprintf("\nDone! Plot saved to:\n  %s\n", out_file))

cat("\n--- Region summary ---\n")
cat(sprintf("Region: Chr%d: %.2f - %.2f Mb\n", chr_num, xlim[1], xlim[2]))
cat(sprintf("Gene midpoint: %.4f Mb\n", gene_mid / 1e6))
cat(sprintf("Lead GWAS SNP: %.4f Mb (p = %.3g)\n", lead_snp / 1e6, lead_p))
cat(sprintf("Max -log10(P) in window: %.2f\n", max(gwas_sub$logP, na.rm = TRUE)))
cat(sprintf("Max FST in window: %.4f\n", max(fst_sub$WEIGHTED_FST, na.rm = TRUE)))
cat(sprintf("FST threshold (99%%): %.4f\n", fst_thresh))
