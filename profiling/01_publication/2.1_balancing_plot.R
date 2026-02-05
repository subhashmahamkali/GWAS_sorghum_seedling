#install.packages("ggrepel")
library(data.table)
library(R.utils)
library(readxl)
library(writexl)
library(dplyr)

#d=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/all_chr_merged_with_chr.B2.txt.gz", header=T,data.table=F)
#d=fread("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/wild_merged_chr.B2.txt", header=T,data.table=F)
d=fread("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/wild_merged_chr.B2.txt", header=T,data.table=F)


d=d[,c(7,1,3)]
#ch=fread("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/sorg.chr_length_V5.txt", header=T,data.table=F)
#ch=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt", header=T,data.table=F)
ch=fread("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/sorg.chr_length_V5.txt", header=T,data.table=F)

d[,1]  <- as.numeric(d[,1])   # chromosome column
d[,2]  <- as.numeric(d[,2])   # position column
d[,3]  <- as.numeric(d[,3])
ch[,4] <- as.numeric(ch[,4])

#gene = read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/W_L_I_5kb_description.xlsx", sheet = 2)
#gene = read_excel("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb_description.xlsx", sheet = 2)
gene = read_excel("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb_description.xlsx", sheet = 2)


gene$pos <- as.integer(gene$pos)
colnames(gene)[1] <- "chr" 
gene$chr <- as.numeric(gene$chr)

for (k in 1:10) {
  sub <- subset(gene, gene$chr == k)
  sub$pos <- sub$pos + ch[k, 4]
  gene[gene$chr == k, "pos"] <- sub$pos
}
gene$pos <- gene$pos / 1e6
#write_xlsx(gene, "W_L_I_5kb_description_processed.xlsx")
red_points3 <- gene[gene$group == "imp", ]  
red_points2 <- gene[gene$group == "lan", ]
red_points1 <- gene[gene$group == "wild", ]

dd=NULL
for (k in 1:10){sub=subset(d,d[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
dd=rbind(dd,sub)}
dd$physPos <- dd[,2] / 1e6 
midpoints <- sapply(1:10, function(k) {
  mean(range(dd[dd[,1] == k, 2], na.rm=TRUE)) / 1e6
})
dd_filt <- dd[dd[,3] >= 0, ]
thr <- quantile(dd_filt[,3], 0.99, na.rm = TRUE)
thr
#831.5203

#la=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)
#la=fread("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)
la=fread("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/Sorghum_landrace.AGPV5.B2_stat.txt", header=T,data.table=F)


la=la[,c(1,2,4)] 
#thr <- quantile(la[,3], 0.99, na.rm = TRUE)
#thr
#1111.447
lan=NULL
for (k in 1:10){sub=subset(la,la[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
lan=rbind(lan,sub)}
lan$physPos <- lan[,2] / 1e6 
midpoins_l <- sapply(1:10, function(k) {
  mean(range(lan[lan[,1] == k, 2], na.rm=TRUE)) / 1e6
})
la_filt <- lan[lan[,3] >= 0, ]
thr_l <- quantile(la_filt[,3], 0.99, na.rm = TRUE)
thr_l
#1111.447

#im=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/imp_chr_merged_with_chr.B2.txt.gz", header=T,data.table=F)
#im=fread("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/imp_merged_chr.B2.txt", header=T,data.table=F)
im=fread("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/imp_merged_chr.B2.txt", header=T,data.table=F)

im=im[,c(7,1,3)]
im[,3] <- as.numeric(im[,3])
#thr <- quantile(im[,3], 0.99, na.rm = TRUE)
#thr
#1682.501
#ch[,4]  <- as.numeric(ch[,4])
im[,1]  <- as.numeric(im[,1])   # chromosome column
im[,2]  <- as.numeric(im[,2])   # position column
im[,3]  <- as.numeric(im[,3])
iim=NULL
for (k in 1:10){sub=subset(im,im[,1]==k)
sub[,2]=sub[,2]+ch[k,4]
iim=rbind(iim,sub)}
#adjust the genomic positiosn for all chromosomes continuesly 
iim$physPos <- iim[,2] / 1e6 
#converts basepairs to megabase
midpoins_i <- sapply(1:10, function(k) {
  mean(range(iim[iim[,1] == k, 2], na.rm=TRUE)) / 1e6
})
ii_filt <- iim[iim[,3] >= 0, ]
thr_i <- quantile(ii_filt[,3], 0.99, na.rm = TRUE)
thr_i
#1682.746

#c = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/W_L_I_5kb.bed")
#c = fread("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb.bed")
c = fread("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb.bed")

w = c[,c(1:5)]
w$V6 <- (w$V2 + w$V3)/2 
w_plot <- NULL
for (k in 1:10) {
  sub = subset(w, w[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V6 =sub$V6 + ch[k,4]
  w_plot = rbind(w_plot, sub)
}
w_plot$V6 = w_plot$V6/1e6 

l = c[,c(6:10)]
l$V11 <- (l$V7 + l$V8)/2 
l_plot <- NULL
for (k in 1:10) {
  sub = subset(l, l[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V11 = sub$V11 + ch[k,4]
  l_plot = rbind(l_plot, sub)
}
l_plot$V11 = l_plot$V11/1e6 

i = c[,c(11:15)]
i$V16 <- (i$V12 + i$V13)/2 
i_plot <- NULL
for (k in 1:10) {
  sub = subset(i, i[[1]] == k)   # w[[1]] is the first column as a vector
  sub$V16 = sub$V16 + ch[k,4]
  i_plot = rbind(i_plot, sub)
}
i_plot$V16 = i_plot$V16/1e6 


# make tidy copies with consistent names
wild_df <- within(dd_filt,  { chr <- dd_filt[,1]; pos_bp <- dd_filt[,2];  B2 <- dd_filt[,3] })
land_df <- within(la_filt,  { chr <- la_filt[,1]; pos_bp <- la_filt[,2];  B2 <- la_filt[,3] })
imp_df  <- within(ii_filt,  { chr <- ii_filt[,1]; pos_bp <- ii_filt[,2];  B2 <- ii_filt[,3] })

# (optional) keep only the standardized columns
wild_df <- wild_df[, c("chr","pos_bp","physPos","B2")]
land_df <- land_df[, c("chr","pos_bp","physPos","B2")]
imp_df  <- imp_df[,  c("chr","pos_bp","physPos","B2")]

chr_id   <- 9
win_loMb <- 601.68708
win_hiMb <- 10913953/1e6
win_mid  <- 607.14406
#gene = read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/W_L_I_5kb_description.xlsx", sheet = 1)
#gene = read_excel("/work/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb_description.xlsx", sheet = 1)
gene = read_excel("/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling/largedata/W_L_I_5kb_description.xlsx", sheet = 1)

gene$pos <- gene$V17+gene$V18/2
gene$V18 <- as.integer(gene$V18)
colnames(gene)[1] <- "chr" 
gene$chr <- as.numeric(gene$chr)

for (k in 1:10) {
  sub <- subset(gene, gene$chr == k)
  sub$V18 <- sub$V18 + ch[k, 4]
  gene[gene$chr == k, "V18"] <- sub$V18
}
gene$st <- gene$V17 / 1e6


sub_w <- subset(wild_df, chr==chr_id & physPos>=win_loMb & physPos<=win_hiMb)
sub_l <- subset(land_df, chr==chr_id & physPos>=win_loMb & physPos<=win_hiMb)
sub_i <- subset(imp_df,  chr==chr_id & physPos>=win_loMb & physPos<=win_hiMb)

gdsl_rows <- tolower(gene$ANNOTATION)
gdsl_pos_candidates <- unique(gene$pos[grepl("^gdsl", gdsl_rows)])  # Mb

## If there are multiple entries, pick the one closest to the center of the zoom
if (length(gdsl_pos_candidates) > 1) {
  gdsl_pos <- gdsl_pos_candidates[ which.min(abs(gdsl_pos_candidates - win_mid)) ]
} else if (length(gdsl_pos_candidates) == 1) {
  gdsl_pos <- gdsl_pos_candidates
} else {
  gdsl_pos <- NA_real_
}

SHOW_GENE_LABELS <- FALSE  # keep dots, hide gene name text
png("graphs/01_publication/2.bal_pos_taj_pi/balancing_plot.png", height = 7, width = 14, res = 600, units = "in")

tiff(
  "graphs/01_publication/2.bal_pos_taj_pi/balancing_plot.tiff",
  width = 14, height = 7, units = "in", res = 600,
  compression = "lzw"
)

par(mar = c(2, 3.5, 2, 2), mfrow = c(3,1), oma = c(2, 2, 2, 1))
# Shared X range (Mb) used by ALL panels
xlim_all <- c(0, 717)
# Colors
col_filt <- ifelse(dd_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
cat_cols  <- c(ion="#1f78b4", nitrogen="#ff7f00", stressrelated="#33a02c")
vline_col <- adjustcolor("grey20", 0.25)

# Local peak helper (±200 kb)
peak_y <- function(df, center_mb, window_mb = 0.20) {
  phys <- if ("physPos" %in% names(df)) df$physPos else df[,2]  # Mb
  b2   <- df[,3]
  lo <- center_mb - window_mb; hi <- center_mb + window_mb
  idx <- phys >= lo & phys <= hi
  if (!any(idx)) return(NA_real_)
  max(b2[idx], na.rm = TRUE)
}

# Precompute dot heights
red_points1$y_peak <- vapply(red_points1$pos, peak_y, numeric(1), df = dd_filt)
red_points2$y_peak <- vapply(red_points2$pos, peak_y, numeric(1), df = la_filt)
red_points3$y_peak <- vapply(red_points3$pos, peak_y, numeric(1), df = ii_filt)

# Shared x positions for guides
vpos_all <- sort(unique(c(red_points1$pos, red_points2$pos, red_points3$pos)))

## ---------- PANEL 1: WILD ----------
plot(
  if ("physPos" %in% names(dd_filt)) dd_filt$physPos else dd_filt[,2],
  dd_filt[,3],
  col = col_filt, pch = 16, cex = 0.4,
  bty = "l",
  xlim = xlim_all,                         # <- SAME xlim for all panels
  axes = FALSE, xlab = "", ylab = "", font.lab = 1
)

#legend("top",
       #legend = c("Ion signalling", "Nitrogen-related", "Stress-related"),
       #col    = cat_cols[c("ion","nitrogen","stressrelated")],
       #pch    = 16,
       #bty    = "o", box.lwd=1.2, box.col="grey30",
       #bg     = adjustcolor("white", 0.9),
       #cex    = 1.0, horiz=TRUE, inset = c(0, -0.2), xpd = NA)

axis(2, las=2, tck=-.03, labels = FALSE) #cex.axis=1, font.axis=1)
chr_breaks <- seq(xlim_all[1], xlim_all[2], length.out=11)
midpoints  <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
#axis(1, at=midpoints, labels=1:10, cex.axis=1.2, font.axis=1)
axis(1, at=midpoints, labels = FALSE, tck=-.03)

segments(x0=xlim_all[1], x1=xlim_all[2], y0=thr, y1=thr, col="red", lty=2, lwd=2)
#mtext("wild sorghum", side=3, line=0.5, font=1, cex=1, adj=1)
#(n=50)
points(red_points1$pos, red_points1$y_peak,
       col = cat_cols[as.character(red_points1$GROUP)], pch = 16, cex = 1)

if (SHOW_GENE_LABELS) {
  x_leader  <- 5; y0 <- 220; y_step <- 180; max_steps <- 10; lab_cex <- 1.2; lab_font <- 3
  placed <- list()
  hits_any <- function(xmin, xmax, ymin, ymax, boxes) {
    if (!length(boxes)) return(FALSE)
    for (b in boxes) if (!(xmax < b[1] || xmin > b[2] || ymax < b[3] || ymin > b[4])) return(TRUE)
    FALSE
  }
  for (i in 1:nrow(red_points1)) {
    x <- red_points1$pos[i]; y <- red_points1$y_peak[i]
    lbl <- tolower(red_points1$ANNOTATION[i])
    col <- cat_cols[as.character(red_points1$GROUP[i])]
    to_right <- (i %% 2 == 0); x_lab <- if (to_right) x + x_leader else x - x_leader
    step <- 0
    repeat {
      y_lab <- y + y0 + step*y_step
      w <- strwidth(lbl, cex=lab_cex, font=lab_font); h <- strheight(lbl, cex=lab_cex, font=lab_font)
      if (to_right) { xmin <- x_lab; xmax <- x_lab + w } else { xmin <- x_lab - w; xmax <- x_lab }
      ymin <- y_lab - 0.5*h; ymax <- y_lab + 0.5*h
      if (!hits_any(xmin, xmax, ymin, ymax, placed) || step >= max_steps) break
      step <- step + 1
    }
    segments(x0=x, y0=y, x1=x_lab, y1=y_lab - 0.2*h, col="grey40", lwd=1.2)
    text(x=x_lab, y=y_lab, labels=lbl, cex=lab_cex, font=lab_font, col=col,
         pos=if (to_right) 4 else 2, xpd=NA)
    placed[[length(placed)+1]] <- c(xmin, xmax, ymin, ymax)
  }
}

# Capture panel plot rectangles (in device coords) for overlay mapping
plt1 <- par("plt")
# Convert x (Mb) -> device x using panel 1's plot rectangle and shared xlim
map_x_ndc <- function(x, plt, xlim) plt[1] + (x - xlim[1]) / diff(xlim) * (plt[2] - plt[1])
x_ndc <- map_x_ndc(vpos_all, plt1, xlim_all)

## ---------- PANEL 2: LANDRACES ----------
col_la <- ifelse(la_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
plot(
  if ("physPos" %in% names(la_filt)) la_filt$physPos else la_filt[,2],
  la_filt[,3],
  col = col_la, pch = 16, cex = 0.4,
  bty = "l", xlim = xlim_all,
  axes = FALSE, xlab = "", ylab = "", font.lab = 1
)
axis(2, las=2, tck=-.03, labels = FALSE)#cex.axis=1, font.axis=1)
midpoints_l <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
#axis(1, at=midpoints_l, labels=1:10, cex.axis=1.2, font.axis=1)
axis(1, at=midpoints, labels = FALSE, tck=-.03)

segments(x0=xlim_all[1], x1=xlim_all[2], y0=thr_l, y1=thr_l, col="red", lty=2, lwd=2)
#mtext("B2", side=2, line=3, font=1, cex=1.5)
#mtext("landrace", side=3, line=0.5, font=1, cex=1, adj=1)
#(n=107)
plt2 <- par("plt")
## ---------- PANEL 3: IMPROVED ----------
col_ii <- ifelse(ii_filt[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
plot(
  if ("physPos" %in% names(ii_filt)) ii_filt$physPos else ii_filt[,2],
  ii_filt[,3],
  col = col_ii, pch = 16, cex = 0.4,
  bty = "l", xlim = xlim_all,
  axes = FALSE, xlab = "", ylab = "", font.lab = 1
)
axis(2, las=2, tck=-.03, labels = FALSE) #cex.axis=1, font.axis=1)
midpoints_i <- (head(chr_breaks,-1) + tail(chr_breaks,-1)) / 2
#axis(1, at=midpoints_i, labels=1:10, cex.axis=1.2, font.axis=1)
axis(1, at=midpoints, labels = FALSE, tck=-.03)

segments(x0=xlim_all[1], x1=xlim_all[2], y0=thr_i, y1=thr_i, col="red", lty=2, lwd=2)
#mtext("Chromosome", side=1, line=2.8, font=1, cex=1.5)
#mtext("improved", side=3, line=0.5, font=1, cex=1, adj=1)
#(n=129)
plt3 <- par("plt")

## ---------- DEVICE-LEVEL OVERLAY: continuous dotted rails ----------
par(new = TRUE, fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")

# span only the union of the three plot rectangles (no outer margins)
y_lo <- min(plt1[3], plt2[3], plt3[3])
y_hi <- max(plt1[4], plt2[4], plt3[4])

for (x in x_ndc) {
  segments(x, y_lo, x, y_hi, lty = 3, lwd = 0.8, col = vline_col, xpd = NA)
}
dev.off()







######################################################################################## Zoomed-in plots for chr9:601.68708-10.913953 Mb region only
par(mfrow = c(3,1), mar = c(3,4,2,1), oma = c(3,0,0,0))

## helper to draw the vertical line + label if the gene sits inside the window
draw_gdsl <- function(df, ylim_pad = 0.95, label = "gdsl") {
  if (is.finite(gdsl_pos) && gdsl_pos >= win_loMb && gdsl_pos <= win_hiMb) {
    abline(v = gdsl_pos, col = "grey40", lty = 3, lwd = 1.5)
    text(gdsl_pos,
         max(df$B2, na.rm = TRUE) * ylim_pad,
         labels = label, font = 3, cex = 0.9, col = "grey20")
  }
}

## wild
plot(sub_w$physPos, sub_w$B2, type = "l", xlab = "", ylab = "B2",
     main = "Balancing selection: wild population")
abline(h = quantile(wild_df$B2, 0.99, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
draw_gdsl(sub_w)

## landraces
plot(sub_l$physPos, sub_l$B2, type = "l", xlab = "", ylab = "B2",
     main = "Balancing selection: landraces population")
abline(h = quantile(land_df$B2, 0.99, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
draw_gdsl(sub_l)

## improved
plot(sub_i$physPos, sub_i$B2, type = "l", xlab = "Position (Mb)", ylab = "B2",
     main = "Balancing selection: improved population")
abline(h = quantile(imp_df$B2, 0.99, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
draw_gdsl(sub_i)
####################################################################################


vpos_all <- sort(unique(c(red_points1$pos, red_points2$pos, red_points3$pos)))
vline_col <- adjustcolor("grey20", 0.25)  # faint
## === Continuous vertical lines across all three panels ===
# Get gene positions in Mb (same as dots)
x_user <- vpos_all   

# Convert those positions into device coords using the top panel
par(mfg = c(1,1))  
x_ndc <- grconvertX(x_user, from = "user", to = "ndc")

# Figure out vertical extent of top and bottom panels
fig_top <- par("fig")
par(mfg = c(3,1))
fig_bot <- par("fig")

y_top <- fig_top[4]
y_bot <- fig_bot[3]

# Overlay full device
par(new = TRUE, fig = c(0,1,0,1), mar = c(0,0,0,0))
plot.new()
segments(x0 = x_ndc, y0 = y_bot, x1 = x_ndc, y1 = y_top,
         col = vline_col, lty = 3, lwd = 0.8, xpd = NA)
## === end continuous lines ===
