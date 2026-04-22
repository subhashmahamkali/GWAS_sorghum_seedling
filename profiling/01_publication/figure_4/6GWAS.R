library(data.table)
library(grDevices)

# =========================
# helper: merge nearby x-positions
# =========================
merge_close_positions <- function(x, tol = 0.15) {
  x <- sort(x[is.finite(x)])
  if (length(x) <= 1) return(x)
  
  groups <- cumsum(c(TRUE, diff(x) > tol))
  merged <- tapply(x, groups, mean)
  as.numeric(merged)
}

# =========================
# 1. Read files
# =========================
a <- fread(
  "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/all_traits.txt",
  header = TRUE, data.table = FALSE
)
a[,1] <- gsub("_filtered.csv", "", a[,1])
a[,1] <- gsub("_", " ", a[,1])
colnames(a)[1] <- "Trait"

ch <- fread(
  "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg.chr_length_V5.txt",
  header = TRUE, data.table = FALSE
)

d_raw <- fread(
  "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/All_GWAS_sorg.txt",
  header = TRUE, data.table = FALSE
)[, -c(5:6)]

# keep original trait strings
d_raw$Trait_full <- d_raw$Trait

# =========================
# 2. Trait mapping
# =========================
trait_map <- data.frame(
  pattern = c(
    "branchInternodeLength",
    "plantHeight",
    "primaryBranchNumber",
    "tillersPerPlant",
    "daysToFlower",
    "extantLeafNumber",
    "flagLeafLength",
    "flagLeafWidth",
    "leafAngleStandardDeviation",
    "medianLeafAngle",
    "stemDiameterLower",
    "stemDiameterUpper",
    "thirdLeafLength",
    "thirdLeafWidth",
    "panicleGrainWeight",
    "paniclesPerPlot",
    "paniclesPerPlant",
    "rachisDiameterLower",
    "rachisDiameterUpper",
    "rachisLength"
  ),
  Trait_clean = c(
    "Branch internode length",
    "Plant height",
    "Primary branch number",
    "Tillers per plant",
    "Days to flower",
    "Extant leaf number",
    "Flag leaf length",
    "Flag leaf width",
    "Leaf angle standard deviation",
    "Median leaf angle",
    "Stem diameter lower",
    "Stem diameter upper",
    "Third leaf length",
    "Third leaf width",
    "Panicle grain weight",
    "Panicles per plot",
    "Panicles per plant",
    "Rachis diameter lower",
    "Rachis diameter upper",
    "Rachis length"
  ),
  Category = c(
    "Architecture",
    "Architecture",
    "Architecture",
    "Architecture",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Developmental",
    "Panicle",
    "Panicle",
    "Panicle",
    "Panicle",
    "Panicle",
    "Panicle"
  ),
  stringsAsFactors = FALSE
)

d_raw$Trait_clean <- NA
d_raw$Category <- NA

for (i in seq_len(nrow(trait_map))) {
  hit <- grepl(trait_map$pattern[i], d_raw$Trait_full, fixed = TRUE)
  d_raw$Trait_clean[hit] <- trait_map$Trait_clean[i]
  d_raw$Category[hit] <- trait_map$Category[i]
}

# keep only wanted traits
d <- d_raw[!is.na(d_raw$Trait_clean), ]

# =========================
# 3. Nitrogen treatment
# =========================
d$Nitrogen_Treatment <- NA
d$Nitrogen_Treatment[grepl("\\.HN\\.|^HN\\.", d$Trait_full)] <- "HN"
d$Nitrogen_Treatment[grepl("\\.LN\\.|^LN\\.", d$Trait_full)] <- "LN"
d$Nitrogen_Treatment[grepl("\\.NR$|\\.NR\\.|^NR\\.|NR$", d$Trait_full)] <- "NR"

d <- d[!is.na(d$Nitrogen_Treatment), ]

# keep columns
d <- d[, c("CHROM", "POS", "P-value", "Trait_clean", "Nitrogen_Treatment", "Category")]
colnames(d)[4] <- "Trait"

# =========================
# 4. Colors
# =========================
col_vec <- c("slateblue", "darkgreen", "violetred", "gold2", "skyblue")
col_vec <- adjustcolor(col_vec, alpha.f = 0.5)

# =========================
# 5. Chromosome formatting
# =========================
d$CHROM <- gsub("Chr", "", d$CHROM)
d$CHROM <- as.integer(d$CHROM)

dd_list <- list()
for (k in 1:10) {
  sub <- subset(d, CHROM == k)
  if (nrow(sub) > 0) {
    sub$POS <- sub$POS + ch[k, 4]
    dd_list[[length(dd_list) + 1]] <- sub
  }
}
dd <- do.call(rbind, dd_list)
dd$POS_MB <- dd$POS / 1e6

# assign colors
dd$col1 <- NA
dd[dd$Category == "Developmental", "col1"] <- col_vec[2]
dd[dd$Category == "Architecture",  "col1"] <- col_vec[1]
dd[dd$Category == "Panicle",       "col1"] <- col_vec[5]

# =========================
# 6. Read overlap loci
# =========================
d2 <- fread(
  "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/NR_LN_HN_ove.txt",
  header = FALSE, data.table = FALSE
)

d3 <- NULL
for (k in 1:10) {
  sub <- subset(d2, d2[,1] == k)
  sub[,5]  <- sub[,5]  + ch[k,4]
  sub[,13] <- sub[,13] + ch[k,4]
  sub[,21] <- sub[,21] + ch[k,4]
  d3 <- rbind(d3, sub)
}

d3[,5]  <- d3[,5]  / 1e6
d3[,13] <- d3[,13] / 1e6
d3[,21] <- d3[,21] / 1e6

# =========================
# 7. One representative x-position per common locus
#    then merge nearby loci to avoid many duplicate dotted lines
# =========================
overlap_x <- rowMeans(cbind(d3[,5], d3[,13], d3[,21]), na.rm = TRUE)
overlap_x <- overlap_x[is.finite(overlap_x)]
overlap_x <- merge_close_positions(overlap_x, tol = 0.15)   # try 0.25 if still crowded

# =========================
# 8. Tick positions
# =========================
y_ticks <- c(2, 4, 6, 8, 10, 12)
x_ticks <- (ch[,3] - ch[,2] / 2) / 1e6

# =========================
# 9. Plot settings
# =========================
xlim_range <- c(0, 720)
ylim_range <- c(2, 12)
sig_line_y <- 5

draw_vertical_overlap_lines <- function(xpos, ylim_range,
                                        col = adjustcolor("black", alpha.f = 0.75),
                                        lty = 3, lwd = 1.0) {
  segments(
    x0 = xpos,
    y0 = ylim_range[1],
    x1 = xpos,
    y1 = ylim_range[2],
    col = col,
    lty = lty,
    lwd = lwd
  )
}

# =========================
# 10. TIFF output
# =========================
tiff(
  "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/up_fig/GWAS_HN_LN_NR_selected_traits_clean_no_text.tiff",
  width = 14,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw"
)

par(
  mfrow = c(3, 1),
  mar   = c(0.35, 2.5, 0.25, 0.5),
  oma   = c(0.45, 0.3, 0.3, 0.3),
  xaxs  = "i",
  yaxs  = "i"
)

# =========================
# 11. HN panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "HN", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1,
  pch = 16,
  cex = 0.45,
  bty = "l",
  xlim = xlim_range,
  ylim = ylim_range,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
segments(
  xlim_range[1], sig_line_y,
  xlim_range[2], sig_line_y,
  col = "red", lty = 2, lwd = 2
)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.3))

# draw dotted lines LAST so they are visible on top of peaks
draw_vertical_overlap_lines(overlap_x, ylim_range)

# =========================
# 12. LN panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "LN", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1,
  pch = 16,
  cex = 0.45,
  bty = "l",
  xlim = xlim_range,
  ylim = ylim_range,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
segments(
  xlim_range[1], sig_line_y,
  xlim_range[2], sig_line_y,
  col = "red", lty = 2, lwd = 2
)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.3))

draw_vertical_overlap_lines(overlap_x, ylim_range)

# =========================
# 13. NR panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "NR", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1,
  pch = 16,
  cex = 0.45,
  bty = "l",
  xlim = xlim_range,
  ylim = ylim_range,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
axis(1, at = x_ticks, labels = FALSE, tck = -0.03)
segments(
  xlim_range[1], sig_line_y,
  xlim_range[2], sig_line_y,
  col = "red", lty = 2, lwd = 2
)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.3))

draw_vertical_overlap_lines(overlap_x, ylim_range)

dev.off()