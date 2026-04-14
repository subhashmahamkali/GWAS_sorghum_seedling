library(data.table)
library(grDevices)

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

# cumulative positions
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
# 7. Tick positions
# =========================
y_ticks <- c(2, 4, 6, 8, 10, 12)
x_ticks <- (ch[,3] - ch[,2] / 2) / 1e6

# =========================
# 8. TIFF output
# =========================
tiff(
  "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/up_fig/GWAS_HN_LN_NR_selected_traits_clean_no_text.tiff",
  width = 14,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw"
)

par(mar = c(1.5, 2.5, 0.4, 0.5), mfrow = c(3, 1), oma = c(0.3, 0.3, 0.3, 0.3))

# =========================
# 9. HN panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "HN", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1, pch = 16, cex = 0.45,
  bty = "l", xlim = c(0, 720), ylim = c(2, 12),
  axes = FALSE, xlab = "", ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
segments(0, 5, 720, 5, col = "red", lty = 2, lwd = 2)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.3))

points(
  d3[,5], -log10(d3[,4]),
  col = adjustcolor("violetred", alpha.f = 0.7),
  pch = 8, cex = 0.8, type = "h", lwd = 1.3
)

# =========================
# 10. LN panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "LN", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1, pch = 16, cex = 0.45,
  bty = "l", xlim = c(0, 720), ylim = c(2, 12),
  axes = FALSE, xlab = "", ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
segments(0, 5, 720, 5, col = "red", lty = 2, lwd = 2)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.6))

points(
  d3[,13], -log10(d3[,12]),
  col = adjustcolor("violetred", alpha.f = 0.7),
  pch = 8, cex = 0.8, type = "h", lwd = 1.3
)

# =========================
# 11. NR panel
# =========================
d1 <- dd[dd$Nitrogen_Treatment == "NR", ]

plot(
  d1$POS_MB, -log10(d1$`P-value`),
  col = d1$col1, pch = 16, cex = 0.45,
  bty = "l", xlim = c(0, 720), ylim = c(2, 12),
  axes = FALSE, xlab = "", ylab = ""
)

axis(2, at = y_ticks, labels = FALSE, las = 2, tck = -0.03)
axis(1, at = x_ticks, labels = FALSE, tck = -0.03)
segments(0, 5, 720, 5, col = "red", lty = 2, lwd = 2)
abline(v = ch[,4] / 1e6, col = adjustcolor("gray", alpha.f = 0.6))

points(
  d3[,21], -log10(d3[,20]),
  col = adjustcolor("violetred", alpha.f = 0.7),
  pch = 8, cex = 0.8, type = "h", lwd = 1.3
)

dev.off()