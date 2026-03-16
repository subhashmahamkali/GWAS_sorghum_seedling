#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub('^--file=', '', args[grepl('^--file=', args)])
if (!length(script_path)) {
  script_path <- sys.frames()[[1]]$ofile
}
if (!length(script_path)) {
  script_path <- "workflow/figure2_7track.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
figure_dir <- file.path(repo_root, "largedata", "figure_2")

path_or_stop <- function(...) {
  p <- file.path(...)
  if (!file.exists(p)) stop("Missing file: ", p)
  p
}

# Core inputs
chr_path     <- path_or_stop(repo_root, "largedata", "figure_2", "sorg.chr_length_V5.txt")
bal_la_unique_path  <- path_or_stop(repo_root, "largedata", "figure_2", "Landrace_only_bal_description.xlsx")
bal_imp_unique_path <- path_or_stop(repo_root, "largedata", "figure_2", "Improved_only_bal_description.xlsx")
b2_wild_path <- path_or_stop(repo_root, "largedata", "wild_merged_chr.B2.txt")
b2_land_path <- path_or_stop(repo_root, "largedata", "Sorghum_landrace.AGPV5.B2_stat.txt")
b2_imp_path  <- path_or_stop(repo_root, "largedata", "imp_merged_chr.B2.txt")
fst_wl_path  <- path_or_stop(repo_root, "largedata", "figure_2", "sorg_wild_vs_landrace.windowed.weir.fst")
fst_li_path  <- path_or_stop(repo_root, "largedata", "figure_2", "landrace_vs_improved.windowed.weir.fst")
xp_dom_path  <- path_or_stop(repo_root, "largedata", "figure_2", "domestication_landrace_vs_wild.all_chr.tsv")
xp_bre_path  <- path_or_stop(repo_root, "largedata", "figure_2", "breeding_improved_vs_landrace.all_chr.tsv")
wli_ann_path <- path_or_stop(repo_root, "largedata", "figure_2", "W_L_I_5kb_description.xlsx")
fst_ann_path <- path_or_stop(repo_root, "largedata", "figure_2", "W_I_L_positive_selection_new.xlsx")
xp_ann_path  <- path_or_stop(repo_root, "largedata", "figure_2", "genes_annotattion_manhattan_plot.xlsx")
out_tiff     <- file.path(figure_dir, "combined_7track_selection.tiff")

chr_tab <- fread(chr_path, data.table = FALSE)
chr_offsets <- setNames(chr_tab$cumsum2, chr_tab$Chr)

continuous_pos <- function(chr_vec, pos_vec) {
  chr_vec <- as.integer(chr_vec)
  pos_vec <- as.numeric(pos_vec)
  offs <- chr_offsets[as.character(chr_vec)]
  (pos_vec + offs) / 1e6
}

load_b2 <- function(path, idx_chr, idx_pos, idx_stat) {
  d <- fread(path, data.table = FALSE)
  data.frame(
    chr = as.integer(d[[idx_chr]]),
    pos_mb = continuous_pos(d[[idx_chr]], d[[idx_pos]]),
    value = as.numeric(d[[idx_stat]])
  )
}

tidy_fst <- function(path) {
  d <- fread(path, data.table = FALSE)
  if (!"WEIGHTED_FST" %in% names(d)) {
    stop("Expected WEIGHTED_FST column in ", path)
  }
  pos <- floor((d$BIN_START + d$BIN_END) / 2)
  data.frame(
    chr = as.integer(d$CHROM),
    pos_mb = continuous_pos(d$CHROM, pos),
    value = as.numeric(d$WEIGHTED_FST)
  )
}

tidy_xpclr <- function(path) {
  d <- fread(path, data.table = FALSE)
  stopifnot(all(c("chrom", "start", "stop", "xpclr") %in% names(d)))
  mid <- floor((d$start + d$stop) / 2)
  data.frame(
    chr = as.integer(d$chrom),
    pos_mb = continuous_pos(d$chrom, mid),
    value = as.numeric(d$xpclr)
  )
}

normalize_hex <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  mask <- grepl("^#", x)
  x[mask & nchar(x) == 5] <- paste0(x[mask & nchar(x) == 5], "0")
  x[mask & nchar(x) == 6] <- paste0(x[mask & nchar(x) == 6], "0")
  valid <- grepl("^#[0-9a-fA-F]{6}$", x)
  x[!valid] <- NA
  x
}

load_balancing_annotations <- function(path) {
  raw <- read_excel(path, sheet = "annoated_onbala_plot")
  cols <- list(
    chr = raw$V16,
    start = raw$w,
    end = raw$V18,
    wild = raw$V5,
    land = raw$V10,
    imp = raw$V15,
    color = normalize_hex(raw$`...27`)
  )
  cat_cols <- c(
    nitrogen = "#ff7f00",
    ion = "#1f78b4",
    stressrelated = "#33a02c"
  )
  fallback <- cat_cols[tolower(raw$GROUP)]
  cols$color[is.na(cols$color)] <- fallback[is.na(cols$color)]
  df <- data.frame(
    chr = as.integer(cols$chr),
    pos_mb = continuous_pos(cols$chr, (as.numeric(cols$start) + as.numeric(cols$end)) / 2),
    wild = as.numeric(cols$wild),
    land = as.numeric(cols$land),
    imp = as.numeric(cols$imp),
    color = cols$color,
    stringsAsFactors = FALSE
  )
  df[complete.cases(df[c("chr", "pos_mb")]), ]
}

load_unique_balancing <- function(path, sheet = "Sheet2") {
  raw <- read_excel(path, sheet = sheet)
  color <- NULL
  if ("dot_color" %in% names(raw)) {
    color <- raw$dot_color
    if (all(is.na(color))) color <- NULL
  }
  if (is.null(color) && "Unnamed: 11" %in% names(raw)) {
    color <- raw$`Unnamed: 11`
  }
  if (is.null(color)) {
    color <- rep(NA_character_, nrow(raw))
  }
  color <- normalize_hex(color)
  df <- data.frame(
    chr = as.integer(raw$V5),
    pos_mb = continuous_pos(raw$V5, (as.numeric(raw$V6) + as.numeric(raw$V7)) / 2),
    value = as.numeric(raw$V4),
    color = color,
    stringsAsFactors = FALSE
  )
  df[complete.cases(df[c("chr", "pos_mb")]), ]
}

load_fst_annotations <- function(path) {
  raw <- read_excel(path, sheet = "Sheet4")
  data.frame(
    chr = as.integer(raw$V6),
    pos_mb = continuous_pos(raw$V6, (as.numeric(raw$V7) + as.numeric(raw$V8)) / 2),
    value = as.numeric(raw$V4),
    group = as.character(raw$group),
    label = as.character(raw$annotation),
    stringsAsFactors = FALSE
  )
}

load_xp_annotations <- function(path) {
  raw <- read_excel(path)
  data.frame(
    comparison = as.character(raw$comparison),
    chr = as.integer(raw$gene_chr),
    pos_mb = continuous_pos(raw$gene_chr, (as.numeric(raw$gene_start) + as.numeric(raw$gene_end)) / 2),
    value = as.numeric(raw$max_xpclr),
    stringsAsFactors = FALSE
  )
}

wild_df <- load_b2(b2_wild_path, idx_chr = 7, idx_pos = 1, idx_stat = 3)
land_df <- load_b2(b2_land_path, idx_chr = 1, idx_pos = 2, idx_stat = 4)
imp_df  <- load_b2(b2_imp_path,  idx_chr = 7, idx_pos = 1, idx_stat = 3)

fst_wl_df <- tidy_fst(fst_wl_path)
fst_li_df <- tidy_fst(fst_li_path)
xp_dom_df <- tidy_xpclr(xp_dom_path)
xp_bre_df <- tidy_xpclr(xp_bre_path)

bal_ann <- load_balancing_annotations(wli_ann_path)
la_unique <- load_unique_balancing(bal_la_unique_path)
imp_unique <- load_unique_balancing(bal_imp_unique_path)
fst_ann <- load_fst_annotations(fst_ann_path)
xp_ann  <- load_xp_annotations(xp_ann_path)

xlim_all <- c(0, max(c(wild_df$pos_mb, land_df$pos_mb, imp_df$pos_mb), na.rm = TRUE))

clip_rows <- function(df) df[df$pos_mb >= xlim_all[1] & df$pos_mb <= xlim_all[2], ]
bal_ann <- clip_rows(bal_ann)
la_unique <- clip_rows(la_unique)
imp_unique <- clip_rows(imp_unique)
fst_ann <- clip_rows(fst_ann)
xp_ann  <- clip_rows(xp_ann)

stopifnot(nrow(fst_ann[fst_ann$group == "WILS_LANDRACE", ]) >= 4)
stopifnot(nrow(fst_ann[fst_ann$group == "LANDRACE_IMPROVED", ]) >= 2)
stopifnot(nrow(xp_ann[xp_ann$comparison == "domestication_landrace_vs_wild", ]) >= 6)
stopifnot(nrow(xp_ann[xp_ann$comparison == "breeding_improved_vs_landrace", ]) >= 4)

fst_wl_ann <- fst_ann[fst_ann$group == "WILS_LANDRACE", ]
fst_wl_ann <- fst_wl_ann[order(-fst_wl_ann$value), ][1:4, ]
fst_li_ann <- fst_ann[fst_ann$group == "LANDRACE_IMPROVED", ]
fst_li_ann <- fst_li_ann[order(-fst_li_ann$value), ][1:2, ]

xp_dom_ann <- xp_ann[xp_ann$comparison == "domestication_landrace_vs_wild", ]
xp_dom_ann <- xp_dom_ann[order(-xp_dom_ann$value), ][1:6, ]
xp_bre_ann <- xp_ann[xp_ann$comparison == "breeding_improved_vs_landrace", ]
xp_bre_ann <- xp_bre_ann[order(-xp_bre_ann$value), ][1:4, ]

bal_cols <- bal_ann$color
bal_cols[is.na(bal_cols)] <- "#888888"

col_fn <- function(chr_vec) {
  chr_vec <- as.integer(chr_vec)
  ifelse(chr_vec %% 2 == 1, "#00000040", "#BEBEBE80")
}

plot_track <- function(df, threshold, ann_df = NULL, show_x = FALSE, v_guides = NULL) {
  plot(df$pos_mb, df$value,
       col = col_fn(df$chr), pch = 16, cex = 0.25,
       bty = "l", xlim = xlim_all,
       axes = FALSE, xlab = "", ylab = "", main = "")
  axis(2, las = 2, tck = -0.02, labels = FALSE)
  if (show_x) {
    axis(1, las = 1, tck = -0.02, labels = FALSE)
  }
  abline(h = threshold, col = "red", lty = 2, lwd = 1)
  if (!is.null(ann_df) && nrow(ann_df) > 0) {
    points(ann_df$pos_mb, ann_df$value,
           pch = 21, bg = ann_df$color, col = "white",
           cex = 1.25, lwd = 0.6)
  }
  if (!is.null(v_guides) && length(v_guides)) {
    abline(v = v_guides,
           col = adjustcolor("grey30", 0.35),
           lty = 3, lwd = 0.8)
  }
}

vz <- sort(unique(bal_ann$pos_mb))

dir.create(dirname(out_tiff), showWarnings = FALSE, recursive = TRUE)
tiff(out_tiff, width = 14, height = 16, units = "in", res = 600,
     compression = "lzw", type = "cairo", bg = "white")
par(mfrow = c(7, 1), mar = c(0.6, 3.6, 0.4, 0.6), oma = c(3, 0.5, 0.5, 0.5))

plt_list <- vector("list", 7)

plot_track(wild_df, quantile(wild_df$value, 0.99, na.rm = TRUE),
           ann_df = data.frame(pos_mb = bal_ann$pos_mb, value = bal_ann$wild, color = bal_cols),
           v_guides = vz)
plt_list[[1]] <- par("plt")
plot_track(land_df, quantile(land_df$value, 0.99, na.rm = TRUE),
           ann_df = la_unique,
           v_guides = vz)
plt_list[[2]] <- par("plt")
plot_track(imp_df, quantile(imp_df$value, 0.99, na.rm = TRUE),
           ann_df = imp_unique,
           v_guides = vz)
plt_list[[3]] <- par("plt")
plot_track(fst_wl_df, quantile(fst_wl_df$value, 0.99, na.rm = TRUE),
           ann_df = transform(fst_wl_ann, color = "#d62728"))
plt_list[[4]] <- par("plt")
plot_track(xp_dom_df, quantile(xp_dom_df$value, 0.99, na.rm = TRUE),
           ann_df = transform(xp_dom_ann, color = "#d62728"))
plt_list[[5]] <- par("plt")
plot_track(fst_li_df, quantile(fst_li_df$value, 0.99, na.rm = TRUE),
           ann_df = transform(fst_li_ann, color = "#d62728"))
plt_list[[6]] <- par("plt")
plot_track(xp_bre_df, quantile(xp_bre_df$value, 0.99, na.rm = TRUE),
           ann_df = transform(xp_bre_ann, color = "#d62728"), show_x = TRUE)
plt_list[[7]] <- par("plt")

dev.off()
cat("Saved 7-track TIFF to", out_tiff, "\n")
