library(data.table)
library(readxl)

setwd("~/Documents/gwas_sap")

# ══════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════

# B2 data (local)
b2_wild_path <- "/Users/subhashmahamkali/Downloads/1.miscellaneous/all_chr_merged_with_chr.B2.txt.gz"
b2_land_path <- "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/Sorghum_landrace.AGPV5.B2_stat.txt"
b2_imp_path  <- "/Users/subhashmahamkali/Downloads/1.miscellaneous/imp_chr_merged_with_chr.B2.txt.gz"

# FST, XPCLR, chromosome lengths (local)
chr_path      <- "data/2.2_positive_selec/sorg.chr_length_V5.txt"
fst_wl_path   <- "data/2.2_positive_selec/sorg_wild_vs_landrace.windowed.weir.fst"
fst_li_path   <- "data/2.2_positive_selec/landrace_vs_improved.windowed.weir.fst"
xpclr_dom_path <- "data/xpclr/domestication_landrace_vs_wild/domestication_landrace_vs_wild.all_chr.tsv"
xpclr_bre_path <- "data/xpclr/breeding_improved_vs_landrace/breeding_improved_vs_landrace.all_chr.tsv"

# Annotation files (local)
ann_wli_path  <- "data/2.0_balancing_selection/01_b2_candidates/W_L_I_5kb_description.xlsx"
ann_land_path <- "data/2.0_balancing_selection/01_b2_candidates/Landrace_only_bal_description.xlsx"
ann_imp_path  <- "data/2.0_balancing_selection/01_b2_candidates/Improved_only_bal_description.xlsx"
ann_fst_path  <- "data/2.2_positive_selec/W_I_L_positive_selection_new.xlsx"
ann_xpclr_path <- "data/xpclr/genes_annotattion_manhattan_plot.xlsx"

# Output
out_path <- "graphs/01_publication/combined_7track_selection.tiff"

# ══════════════════════════════════════════════════════════════════
# CHROMOSOME LENGTHS
# ══════════════════════════════════════════════════════════════════

ch <- fread(chr_path, header = TRUE, data.table = FALSE)
# columns: Chr, Pos (chr length), cumsum1 (cumulative start), cumsum2 (cumulative offset for plotting)

# ══════════════════════════════════════════════════════════════════
# HELPER: add cumulative offset and convert to Mb
# ══════════════════════════════════════════════════════════════════

add_offset_mb <- function(df, chr_col = 1, pos_col = 2) {
  result <- NULL
  for (k in 1:10) {
    sub <- df[df[[chr_col]] == k, ]
    if (nrow(sub) == 0) next
    sub[[pos_col]] <- sub[[pos_col]] + ch[k, 4]
    result <- rbind(result, sub)
  }
  result[[pos_col]] <- result[[pos_col]] / 1e6
  result
}

# ══════════════════════════════════════════════════════════════════
# LOAD B2 DATA
# ══════════════════════════════════════════════════════════════════

# Wild: cols 7 (chr), 1 (pos), 3 (B2)
wild_raw <- fread(b2_wild_path, header = TRUE, data.table = FALSE)
wild_raw <- wild_raw[, c(7, 1, 3)]
wild_raw[, 1] <- as.numeric(wild_raw[, 1])
wild_raw[, 2] <- as.numeric(wild_raw[, 2])
wild_raw[, 3] <- as.numeric(wild_raw[, 3])
names(wild_raw) <- c("chr", "pos", "B2")
dd <- NULL
for (k in 1:10) {
  sub <- wild_raw[wild_raw$chr == k, ]
  sub$pos <- sub$pos + ch[k, 4]
  dd <- rbind(dd, sub)
}
dd$pos_mb <- dd$pos / 1e6
dd_filt <- dd[!is.na(dd$B2) & dd$B2 >= 0, ]
thr_w <- quantile(dd_filt$B2, 0.99, na.rm = TRUE)

# Landrace: cols 1 (chr), 2 (pos), 4 (B2)
la_raw <- fread(b2_land_path, header = TRUE, data.table = FALSE)
la_raw <- la_raw[, c(1, 2, 4)]
names(la_raw) <- c("chr", "pos", "B2")
lan <- NULL
for (k in 1:10) {
  sub <- la_raw[la_raw$chr == k, ]
  sub$pos <- sub$pos + ch[k, 4]
  lan <- rbind(lan, sub)
}
lan$pos_mb <- lan$pos / 1e6
la_filt <- lan[!is.na(lan$B2) & lan$B2 >= 0, ]
thr_l <- quantile(la_filt$B2, 0.99, na.rm = TRUE)

# Improved: cols 7 (chr), 1 (pos), 3 (B2)
im_raw <- fread(b2_imp_path, header = TRUE, data.table = FALSE)
im_raw <- im_raw[, c(7, 1, 3)]
im_raw[, 1] <- as.numeric(im_raw[, 1])
im_raw[, 2] <- as.numeric(im_raw[, 2])
im_raw[, 3] <- as.numeric(im_raw[, 3])
names(im_raw) <- c("chr", "pos", "B2")
iim <- NULL
for (k in 1:10) {
  sub <- im_raw[im_raw$chr == k, ]
  sub$pos <- sub$pos + ch[k, 4]
  iim <- rbind(iim, sub)
}
iim$pos_mb <- iim$pos / 1e6
ii_filt <- iim[!is.na(iim$B2) & iim$B2 >= 0, ]
thr_i <- quantile(ii_filt$B2, 0.99, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════
# LOAD FST DATA
# ══════════════════════════════════════════════════════════════════

prep_fst <- function(path) {
  d <- fread(path, header = TRUE, data.table = FALSE)[, -c(4, 6)]
  d$POS <- (d[, 2] + d[, 3]) / 2
  d$POS <- ifelse(d$POS %% 1 == 0.5, floor(d$POS), d$POS)
  d <- d[, c(1, 3, 2)]               # chr, FST, midpoint
  d <- d[order(d[, 1], d[, 2]), ]
  d <- d[!is.na(d[, 2]), ]
  dp <- NULL
  for (k in 1:10) {
    sub <- d[d[, 1] == k, ]
    sub[, 2] <- sub[, 2] + ch[k, 4]
    dp <- rbind(dp, sub)
  }
  dp[, 2] <- dp[, 2] / 1e6
  dp
}

dp1 <- prep_fst(fst_wl_path)    # Wild vs Landrace (domestication)
dp2 <- prep_fst(fst_li_path)    # Landrace vs Improved (breeding)
thr_fst1 <- quantile(dp1[, 3], 0.99, na.rm = TRUE)
thr_fst2 <- quantile(dp2[, 3], 0.99, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════
# LOAD XPCLR DATA
# ══════════════════════════════════════════════════════════════════

prep_xpclr <- function(path) {
  d <- fread(path, data.table = FALSE)
  d$chrom <- as.integer(d$chrom)
  d$mid <- floor((d$start + d$stop) / 2)
  dp <- NULL
  for (k in 1:10) {
    sub <- d[d$chrom == k, ]
    sub$mid <- sub$mid + ch[k, 4]
    dp <- rbind(dp, sub)
  }
  dp$pos_mb <- dp$mid / 1e6
  dp[is.finite(dp$xpclr) & is.finite(dp$pos_mb), ]
}

xp_dom <- prep_xpclr(xpclr_dom_path)
xp_bre <- prep_xpclr(xpclr_bre_path)
thr_xp_dom <- quantile(xp_dom$xpclr, 0.99, na.rm = TRUE)
thr_xp_bre <- quantile(xp_bre$xpclr, 0.99, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════
# LOAD ANNOTATION DATA
# ══════════════════════════════════════════════════════════════════

# Helper: fix truncated hex colors from Excel
fix_color <- function(x) {
  x <- as.character(x)
  # pad to 7 chars (#rrggbb) if truncated
  x <- ifelse(nchar(x) == 6 & startsWith(x, "#"), paste0(x, "0"), x)
  x <- ifelse(nchar(x) == 5 & startsWith(x, "#"), paste0(x, "b4"), x)  # #1f78b -> #1f78b4
  x
}

# Helper: gene midpoint -> cumulative Mb
gene_to_mb <- function(chr_vec, start_vec, end_vec) {
  mid_bp <- (as.numeric(start_vec) + as.numeric(end_vec)) / 2
  for (k in 1:10) {
    idx <- as.numeric(chr_vec) == k
    mid_bp[idx] <- mid_bp[idx] + ch[k, 4]
  }
  mid_bp / 1e6
}

# ─── 1. Shared W_L_I loci (dots on all 3 B2 panels + vertical lines) ───
wli_raw <- read_excel(ann_wli_path, sheet = "annoated_onbala_plot")
wli_ann <- data.frame(
  gene_chr  = as.numeric(wli_raw[[19]]),
  gene_start = as.numeric(wli_raw[[20]]),
  gene_end   = as.numeric(wli_raw[[21]]),
  wild_B2   = as.numeric(wli_raw[[5]]),
  land_B2   = as.numeric(wli_raw[[11]]),
  imp_B2    = as.numeric(wli_raw[[17]]),
  dot_color = fix_color(wli_raw[[27]]),
  stringsAsFactors = FALSE
)
wli_ann <- wli_ann[!is.na(wli_ann$gene_chr), ]
wli_ann$pos_mb <- gene_to_mb(wli_ann$gene_chr, wli_ann$gene_start, wli_ann$gene_end)

# ─── 2. Landrace-only dots (landrace panel only, no lines) ───
la_only_raw <- read_excel(ann_land_path, sheet = 2)
la_ann <- data.frame(
  gene_chr   = as.numeric(la_only_raw[[5]]),
  gene_start = as.numeric(la_only_raw[[6]]),
  gene_end   = as.numeric(la_only_raw[[7]]),
  land_B2    = as.numeric(la_only_raw[[4]]),
  dot_color  = fix_color(la_only_raw[[ncol(la_only_raw)]]),
  stringsAsFactors = FALSE
)
la_ann <- la_ann[!is.na(la_ann$gene_chr) & !is.na(la_ann$dot_color) & la_ann$dot_color != "NA", ]
la_ann$pos_mb <- gene_to_mb(la_ann$gene_chr, la_ann$gene_start, la_ann$gene_end)

# ─── 3. Improved-only dots (improved panel only, no lines) ───
imp_only_raw <- read_excel(ann_imp_path, sheet = 2)
imp_ann <- data.frame(
  gene_chr   = as.numeric(imp_only_raw[[5]]),
  gene_start = as.numeric(imp_only_raw[[6]]),
  gene_end   = as.numeric(imp_only_raw[[7]]),
  imp_B2     = as.numeric(imp_only_raw[[4]]),
  dot_color  = fix_color(imp_only_raw[[11]]),
  stringsAsFactors = FALSE
)
imp_ann <- imp_ann[!is.na(imp_ann$gene_chr) & !is.na(imp_ann$dot_color), ]
imp_ann$pos_mb <- gene_to_mb(imp_ann$gene_chr, imp_ann$gene_start, imp_ann$gene_end)

# ─── 4. FST annotation dots ───
fst_raw <- read_excel(ann_fst_path, sheet = "Sheet4")
fst_ann <- data.frame(
  gene_chr   = as.numeric(fst_raw[[6]]),
  gene_start = as.numeric(fst_raw[[7]]),
  gene_end   = as.numeric(fst_raw[[8]]),
  fst_val    = as.numeric(fst_raw[[4]]),
  group      = as.character(fst_raw[[12]]),
  stringsAsFactors = FALSE
)
fst_ann <- fst_ann[!is.na(fst_ann$gene_chr), ]
fst_ann$pos_mb <- gene_to_mb(fst_ann$gene_chr, fst_ann$gene_start, fst_ann$gene_end)
fst_wl_ann <- fst_ann[fst_ann$group == "WILS_LANDRACE", ]
fst_li_ann <- fst_ann[fst_ann$group == "LANDRACE_IMPROVED", ]

# ─── 5. XPCLR annotation dots ───
xp_raw <- read_excel(ann_xpclr_path, sheet = 1)
xp_raw <- xp_raw[!is.na(xp_raw[[1]]), ]   # drop empty trailing row
xp_ann <- data.frame(
  comparison = as.character(xp_raw[[1]]),
  gene_chr   = as.numeric(xp_raw[[7]]),
  gene_start = as.numeric(xp_raw[[8]]),
  gene_end   = as.numeric(xp_raw[[9]]),
  xpclr_val  = as.numeric(xp_raw[[5]]),
  stringsAsFactors = FALSE
)
xp_ann <- xp_ann[!is.na(xp_ann$gene_chr), ]
xp_ann$pos_mb <- gene_to_mb(xp_ann$gene_chr, xp_ann$gene_start, xp_ann$gene_end)
xp_dom_ann <- xp_ann[xp_ann$comparison == "domestication_landrace_vs_wild", ]
xp_bre_ann <- xp_ann[xp_ann$comparison == "breeding_improved_vs_landrace", ]

# ══════════════════════════════════════════════════════════════════
# PLOT
# ══════════════════════════════════════════════════════════════════

xlim_all <- c(0, 717)
col_fn   <- function(chr_vec) ifelse(as.integer(chr_vec) %% 2 == 1, "#00000066", "#BEBEBE99")

tiff(out_path, width = 14, height = 16, units = "in", res = 600,
     compression = "lzw", type = "cairo", bg = "white")
par(mfrow = c(7, 1),
    mar   = c(0.25, 3, 0.25, 0.5),
    oma   = c(3, 0.5, 0.5, 0.5))

plts <- vector("list", 7)   # store plt coords for vertical lines

# ─── Panel drawing helper ────────────────────────────────────────
draw_panel <- function(x, y, chr_vec, threshold,
                       ann_x = NULL, ann_y = NULL, ann_col = NULL,
                       is_bottom = FALSE) {
  plot(x, y,
       col  = col_fn(chr_vec),
       pch  = 16, cex = 0.3,
       bty  = "l",
       xlim = xlim_all,
       axes = FALSE, xlab = "", ylab = "")

  axis(2, las = 2, tck = -0.025, labels = FALSE)
  if (is_bottom) axis(1, tck = -0.025, labels = FALSE)

  segments(x0 = xlim_all[1], x1 = xlim_all[2],
           y0 = threshold,   y1 = threshold,
           col = "red", lty = 2, lwd = 1.5)

  if (!is.null(ann_x) && length(ann_x) > 0) {
    points(ann_x, ann_y,
           pch = 21, bg = ann_col, col = "white",
           cex = 1.4, lwd = 0.6)
  }
}

# ─── Panel 1: B2 Wild ────────────────────────────────────────────
draw_panel(dd_filt$pos_mb, dd_filt$B2, dd_filt$chr, thr_w,
           ann_x   = wli_ann$pos_mb,
           ann_y   = wli_ann$wild_B2,
           ann_col = wli_ann$dot_color)
plts[[1]] <- par("plt")

# ─── Panel 2: B2 Landrace ────────────────────────────────────────
draw_panel(la_filt$pos_mb, la_filt$B2, la_filt$chr, thr_l)
# shared loci dots
points(wli_ann$pos_mb, wli_ann$land_B2,
       pch = 21, bg = wli_ann$dot_color, col = "white", cex = 1.4, lwd = 0.6)
# landrace-only dots
points(la_ann$pos_mb, la_ann$land_B2,
       pch = 21, bg = la_ann$dot_color, col = "white", cex = 1.4, lwd = 0.6)
plts[[2]] <- par("plt")

# ─── Panel 3: B2 Improved ────────────────────────────────────────
draw_panel(ii_filt$pos_mb, ii_filt$B2, ii_filt$chr, thr_i)
# shared loci dots
points(wli_ann$pos_mb, wli_ann$imp_B2,
       pch = 21, bg = wli_ann$dot_color, col = "white", cex = 1.4, lwd = 0.6)
# improved-only dots
points(imp_ann$pos_mb, imp_ann$imp_B2,
       pch = 21, bg = imp_ann$dot_color, col = "white", cex = 1.4, lwd = 0.6)
plts[[3]] <- par("plt")

# ─── Panel 4: FST Wild vs Landrace (domestication) ───────────────
draw_panel(dp1[, 2], dp1[, 3], dp1[, 1], thr_fst1,
           ann_x   = fst_wl_ann$pos_mb,
           ann_y   = fst_wl_ann$fst_val,
           ann_col = rep("red", nrow(fst_wl_ann)))
plts[[4]] <- par("plt")

# ─── Panel 5: XPCLR Domestication ───────────────────────────────
draw_panel(xp_dom$pos_mb, xp_dom$xpclr, xp_dom$chrom, thr_xp_dom,
           ann_x   = xp_dom_ann$pos_mb,
           ann_y   = xp_dom_ann$xpclr_val,
           ann_col = rep("red", nrow(xp_dom_ann)))
plts[[5]] <- par("plt")

# ─── Panel 6: FST Landrace vs Improved (breeding) ────────────────
draw_panel(dp2[, 2], dp2[, 3], dp2[, 1], thr_fst2,
           ann_x   = fst_li_ann$pos_mb,
           ann_y   = fst_li_ann$fst_val,
           ann_col = rep("red", nrow(fst_li_ann)))
plts[[6]] <- par("plt")

# ─── Panel 7: XPCLR Breeding (bottom — x-axis ticks) ────────────
draw_panel(xp_bre$pos_mb, xp_bre$xpclr, xp_bre$chrom, thr_xp_bre,
           ann_x     = xp_bre_ann$pos_mb,
           ann_y     = xp_bre_ann$xpclr_val,
           ann_col   = rep("red", nrow(xp_bre_ann)),
           is_bottom = TRUE)
plts[[7]] <- par("plt")

# ══════════════════════════════════════════════════════════════════
# DEVICE-LEVEL VERTICAL DOTTED LINES
# (shared W_L_I loci connecting all 3 B2 panels only)
# ══════════════════════════════════════════════════════════════════

# Get unique gene positions (one line per gene in the shared set)
vpos_mb <- sort(unique(wli_ann$pos_mb))

# Convert from panel-1 user coords to NDC
par(mfg = c(1, 1))
x_ndc <- grconvertX(vpos_mb, from = "user", to = "ndc")

# Vertical span: top of panel 1 → bottom of panel 3
y_top <- plts[[1]][4]
y_bot <- plts[[3]][3]

# Draw across device
par(new = TRUE, fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
segments(x0 = x_ndc, y0 = y_bot, x1 = x_ndc, y1 = y_top,
         lty = 3, lwd = 0.8,
         col = adjustcolor("grey30", 0.35),
         xpd = NA)

dev.off()
message("Done: ", out_path)
