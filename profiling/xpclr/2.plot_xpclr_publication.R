library(data.table)
library(readxl)

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
annot_file <- file.path(repo_root, "data/xpclr/genes_annotattion_manhattan_plot.xlsx")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ch <- fread(chr_file)

annot <- read_excel(annot_file, sheet = 1) |> as.data.frame()
annot$region_chr <- as.integer(annot$region_chr)
annot$region_start <- as.integer(annot$region_start)
annot$region_end <- as.integer(annot$region_end)
annot$max_xpclr <- as.numeric(annot$max_xpclr)
annot$region_mid <- floor((annot$region_start + annot$region_end) / 2)
annot$pos_cum <- annot$region_mid + ch$cumsum2[match(annot$region_chr, ch$Chr)]
annot$pos_mb <- annot$pos_cum / 1e6

prepare_xpclr <- function(path) {
  d <- fread(path)
  stopifnot(all(c("chrom", "start", "stop", "xpclr") %in% names(d)))

  d[, chrom := as.integer(chrom)]
  d[, pos := (start + stop) / 2]
  d[, pos := ifelse(pos %% 1 == 0.5, floor(pos), pos)]
  setorder(d, chrom, pos)

  dp <- NULL
  for (k in 1:10) {
    sub <- d[chrom == k]
    sub[, pos_cum := pos + ch[k, cumsum2]]
    dp <- rbind(dp, sub, fill = TRUE)
  }

  dp[, pos_mb := pos_cum / 1e6]
  dp[, col := ifelse(chrom %% 2 == 1, "#00000066", "#BEBEBE99")]

  list(
    data = dp,
    thr = as.numeric(quantile(dp$xpclr, 0.99, na.rm = TRUE))
  )
}

dom <- prepare_xpclr(dom_file)
breed <- prepare_xpclr(breed_file)

x_ticks <- (ch$cumsum2 + ch$Pos / 2) / 1e6
x_max <- 717

plot_panel <- function(dp, thr, panel_title, label_key, show_x_axis = FALSE) {
  plot(
    dp$pos_mb, dp$xpclr,
    col = dp$col,
    pch = 16,
    cex = 0.4,
    bty = "l",
    xlim = c(0, x_max),
    axes = FALSE,
    cex.lab = 0.6,
    xlab = "",
    ylab = "",
    font.lab = 2
  )

  axis(2, las = 2, tck = -0.03, cex.axis = 1.1, font.axis = 1)
  axis(
    1,
    at = x_ticks,
    tck = -0.03,
    cex.axis = 1.2,
    labels = if (show_x_axis) 1:10 else FALSE
  )
  segments(x0 = 0, x1 = x_max, y0 = thr, y1 = thr, col = "red", lty = 2, lwd = 2)

  mark <- subset(annot, comparison == label_key)
  if (nrow(mark) > 0) {
    segments(
      x0 = mark$pos_mb,
      x1 = mark$pos_mb,
      y0 = 0,
      y1 = mark$max_xpclr,
      col = "#FF69B466",
      lwd = 0.8
    )
  }

  mtext(panel_title, side = 3, line = 0.2, adj = 0.98, font = 2, cex = 1)
  mtext("XP-CLR", side = 2, line = 2.8, font = 2, cex = 1)
  if (show_x_axis) {
    mtext("Chromosome", side = 1, line = 2.2, font = 2, cex = 1.1)
  }
}

tiff(
  file.path(out_dir, "xpclr_positive_selection.tiff"),
  height = 7,
  width = 14,
  res = 600,
  units = "in",
  compression = "lzw",
  type = "cairo",
  bg = "white"
)

par(mar = c(1.5, 4.5, 2.2, 1.5), mfrow = c(2, 1), oma = c(1.5, 0.5, 0.5, 0.5))
plot_panel(dom$data, dom$thr, "Landrace vs Wild", "domestication_landrace_vs_wild", show_x_axis = FALSE)
plot_panel(breed$data, breed$thr, "Improved vs Landrace", "breeding_improved_vs_landrace", show_x_axis = TRUE)
dev.off()

pdf(
  file.path(out_dir, "xpclr_positive_selection.pdf"),
  width = 14,
  height = 7
)

par(mar = c(1.5, 4.5, 2.2, 1.5), mfrow = c(2, 1), oma = c(1.5, 0.5, 0.5, 0.5))
plot_panel(dom$data, dom$thr, "Landrace vs Wild", "domestication_landrace_vs_wild", show_x_axis = FALSE)
plot_panel(breed$data, breed$thr, "Improved vs Landrace", "breeding_improved_vs_landrace", show_x_axis = TRUE)
dev.off()

thresholds <- data.table(
  comparison = c("domestication_landrace_vs_wild", "breeding_improved_vs_landrace"),
  top1pct_threshold = c(dom$thr, breed$thr)
)

fwrite(
  thresholds,
  file.path(out_dir, "xpclr_top1pct_thresholds.txt"),
  sep = "\t"
)

annot_minimal <- subset(
  annot,
  (comparison == "domestication_landrace_vs_wild" &
     region_chr %in% c(1L, 4L, 6L) &
     paste(region_start, region_end, sep = "-") %in% c(
       "4410001-4540000",
       "9555001-9640000",
       "63965001-64085000",
       "70035001-70315000",
       "66615001-66695000",
       "48970001-49155000"
     )) |
    (comparison == "breeding_improved_vs_landrace" &
       region_chr %in% c(1L, 2L, 4L, 7L) &
       paste(region_start, region_end, sep = "-") %in% c(
         "9260001-9335000",
         "76640001-76770000",
         "320001-345000",
         "59645001-59705000"
       ))
)

plot_panel_minimal <- function(dp, label_key) {
  plot(
    dp$pos_mb, dp$xpclr,
    col = dp$col,
    pch = 16,
    cex = 0.4,
    bty = "l",
    xlim = c(0, x_max),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  axis(1, labels = FALSE, tck = -0.02)
  axis(2, labels = FALSE, tck = -0.02, las = 2)

  mark <- subset(annot_minimal, comparison == label_key)
  if (nrow(mark) > 0) {
    segments(
      x0 = mark$pos_mb,
      x1 = mark$pos_mb,
      y0 = 0,
      y1 = mark$max_xpclr,
      col = "#FF69B466",
      lwd = 0.8
    )
  }
}

tiff(
  file.path(out_dir, "xpclr_positive_selection_minimal_annotations.tiff"),
  height = 7,
  width = 14,
  res = 600,
  units = "in",
  compression = "lzw",
  type = "cairo",
  bg = "white"
)

par(mar = c(0.4, 0.4, 0.4, 0.4), mfrow = c(2, 1), oma = c(0, 0, 0, 0))
plot_panel_minimal(dom$data, "domestication_landrace_vs_wild")
plot_panel_minimal(breed$data, "breeding_improved_vs_landrace")
dev.off()
