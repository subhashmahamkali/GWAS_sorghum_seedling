suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# =========================
# Interactive paths
# =========================
repo_root  <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling"
figure_dir <- file.path(repo_root, "largedata", "figure_2")

path_or_stop <- function(...) {
  p <- file.path(...)
  if (!file.exists(p)) stop("Missing file: ", p)
  p
}

# =========================
# Input files
# =========================
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

out_tiff <- file.path(figure_dir, "combined_7track_selection_interactive_final_no_transparency_compact_fixedticks.tiff")

# =========================
# Helpers
# =========================
chr_tab <- fread(chr_path, data.table = FALSE)

# File columns are: Chr, Pos, cumsum1, cumsum2
# Pos is chromosome length
chr_offsets <- setNames(chr_tab$cumsum2, chr_tab$Chr)
chr_centers_mb <- (chr_tab$cumsum2 + chr_tab$Pos / 2) / 1e6

continuous_pos <- function(chr_vec, pos_vec) {
  chr_vec <- as.integer(chr_vec)
  pos_vec <- as.numeric(pos_vec)
  offs <- chr_offsets[as.character(chr_vec)]
  (pos_vec + offs) / 1e6
}

normalize_hex <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "FALSE")] <- NA
  x <- toupper(x)
  out <- rep(NA_character_, length(x))
  keep <- grepl("^#[0-9A-F]{6}$", x)
  out[keep] <- x[keep]
  out
}

# No transparency
col_fn <- function(chr_vec) {
  chr_vec <- as.integer(chr_vec)
  ifelse(chr_vec %% 2 == 1, "#000000", "#BEBEBE")
}

clip_rows <- function(df, xlim_all) {
  df[df$pos_mb >= xlim_all[1] & df$pos_mb <= xlim_all[2], , drop = FALSE]
}

tick_store <- list()

# Fixed rounded y-ticks for each panel
fixed_ticks <- list(
  wild_B2 = c(0, 1200, 2400, 3600),
  landrace_B2 = c(0, 2000, 4000, 6000),
  improved_B2 = c(0, 1500, 3000, 4500),
  fst_wild_landrace = c(0.00, 0.25, 0.50, 0.75),
  xpclr_domestication = c(0, 500, 1000, 1500),
  fst_landrace_improved = c(0.00, 0.15, 0.30, 0.45),
  xpclr_breeding = c(0, 300, 600, 900)
)

# =========================
# Readers
# =========================
load_b2 <- function(path, idx_chr, idx_pos, idx_stat) {
  d <- fread(path, select = c(idx_chr, idx_pos, idx_stat), data.table = FALSE)
  data.frame(
    chr    = as.integer(d[[1]]),
    pos_mb = continuous_pos(d[[1]], d[[2]]),
    value  = as.numeric(d[[3]])
  )
}

tidy_fst <- function(path) {
  d <- fread(
    path,
    select = c("CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST"),
    data.table = FALSE
  )
  pos <- floor((d$BIN_START + d$BIN_END) / 2)
  data.frame(
    chr    = as.integer(d$CHROM),
    pos_mb = continuous_pos(d$CHROM, pos),
    value  = as.numeric(d$WEIGHTED_FST)
  )
}

tidy_xpclr <- function(path) {
  d <- fread(
    path,
    select = c("chrom", "start", "stop", "xpclr"),
    data.table = FALSE
  )
  mid <- floor((d$start + d$stop) / 2)
  data.frame(
    chr    = as.integer(d$chrom),
    pos_mb = continuous_pos(d$chrom, mid),
    value  = as.numeric(d$xpclr)
  )
}

load_balancing_annotations <- function(path) {
  raw <- read_excel(path, sheet = "annoated_onbala_plot")
  
  color <- normalize_hex(raw$`...27`)
  
  cat_cols <- c(
    nitrogen      = "#FF7F00",
    ion           = "#1F78B4",
    stressrelated = "#33A02C"
  )
  
  grp <- tolower(as.character(raw$GROUP))
  fallback <- unname(cat_cols[grp])
  color[is.na(color)] <- fallback[is.na(color)]
  
  df <- data.frame(
    chr    = as.integer(raw$V16),
    pos_mb = continuous_pos(raw$V16, (as.numeric(raw$w) + as.numeric(raw$V18)) / 2),
    wild   = as.numeric(raw$V5),
    land   = as.numeric(raw$V10),
    imp    = as.numeric(raw$V15),
    color  = color,
    stringsAsFactors = FALSE
  )
  
  df[complete.cases(df[, c("chr", "pos_mb")]), , drop = FALSE]
}

load_unique_balancing <- function(path, sheet = "Sheet2", fallback_colors = NULL, default_color = "#1F78B4") {
  raw <- read_excel(path, sheet = sheet)
  
  color <- NULL
  
  if ("dot_color" %in% names(raw)) {
    color <- normalize_hex(raw$dot_color)
    if (all(is.na(color))) color <- NULL
  }
  
  if (is.null(color) && "...12" %in% names(raw)) {
    color <- normalize_hex(raw$...12)
    if (all(is.na(color))) color <- NULL
  }
  
  if (is.null(color)) {
    color <- rep(NA_character_, nrow(raw))
  }
  
  if (!is.null(fallback_colors)) {
    if (length(fallback_colors) != nrow(raw)) {
      stop("fallback_colors length must match number of rows in ", path)
    }
    fallback_colors <- toupper(fallback_colors)
    color[is.na(color)] <- fallback_colors[is.na(color)]
  }
  
  color[is.na(color)] <- toupper(default_color)
  
  df <- data.frame(
    chr    = as.integer(raw$V5),
    pos_mb = continuous_pos(raw$V5, (as.numeric(raw$V6) + as.numeric(raw$V7)) / 2),
    value  = as.numeric(raw$V4),
    color  = color,
    stringsAsFactors = FALSE
  )
  
  df[complete.cases(df[, c("chr", "pos_mb", "value")]), , drop = FALSE]
}

load_fst_annotations <- function(path) {
  raw <- read_excel(path, sheet = "Sheet4")
  data.frame(
    chr    = as.integer(raw$V6),
    pos_mb = continuous_pos(raw$V6, (as.numeric(raw$V7) + as.numeric(raw$V8)) / 2),
    value  = as.numeric(raw$V4),
    group  = as.character(raw$group),
    label  = as.character(raw$annotation),
    stringsAsFactors = FALSE
  )
}

load_xp_annotations <- function(path) {
  raw <- read_excel(path)
  data.frame(
    comparison = as.character(raw$comparison),
    chr        = as.integer(raw$gene_chr),
    pos_mb     = continuous_pos(raw$gene_chr, (as.numeric(raw$gene_start) + as.numeric(raw$gene_end)) / 2),
    value      = as.numeric(raw$max_xpclr),
    stringsAsFactors = FALSE
  )
}

# =========================
# Load data
# =========================
wild_df <- load_b2(b2_wild_path, idx_chr = 7, idx_pos = 1, idx_stat = 3)
land_df <- load_b2(b2_land_path, idx_chr = 1, idx_pos = 2, idx_stat = 4)
imp_df  <- load_b2(b2_imp_path,  idx_chr = 7, idx_pos = 1, idx_stat = 3)

fst_wl_df <- tidy_fst(fst_wl_path)
fst_li_df <- tidy_fst(fst_li_path)
xp_dom_df <- tidy_xpclr(xp_dom_path)
xp_bre_df <- tidy_xpclr(xp_bre_path)

bal_ann <- load_balancing_annotations(wli_ann_path)

la_unique <- load_unique_balancing(
  bal_la_unique_path,
  fallback_colors = c("#33A02C", "#33A02C", "#1F78B4")
)

imp_unique <- load_unique_balancing(
  bal_imp_unique_path,
  default_color = "#1F78B4"
)

fst_ann <- load_fst_annotations(fst_ann_path)
xp_ann  <- load_xp_annotations(xp_ann_path)

# =========================
# Clip
# =========================
xlim_all <- c(0, max(c(wild_df$pos_mb, land_df$pos_mb, imp_df$pos_mb), na.rm = TRUE))

bal_ann    <- clip_rows(bal_ann, xlim_all)
la_unique  <- clip_rows(la_unique, xlim_all)
imp_unique <- clip_rows(imp_unique, xlim_all)
fst_ann    <- clip_rows(fst_ann, xlim_all)
xp_ann     <- clip_rows(xp_ann, xlim_all)

# =========================
# Top annotation points
# =========================
fst_wl_ann <- fst_ann[fst_ann$group == "WILS_LANDRACE", , drop = FALSE]
fst_wl_ann <- fst_wl_ann[order(-fst_wl_ann$value), ][1:min(4, nrow(fst_wl_ann)), , drop = FALSE]

fst_li_ann <- fst_ann[fst_ann$group == "LANDRACE_IMPROVED", , drop = FALSE]
fst_li_ann <- fst_li_ann[order(-fst_li_ann$value), ][1:min(2, nrow(fst_li_ann)), , drop = FALSE]

xp_dom_ann <- xp_ann[xp_ann$comparison == "domestication_landrace_vs_wild", , drop = FALSE]
xp_dom_ann <- xp_dom_ann[order(-xp_dom_ann$value), ][1:min(6, nrow(xp_dom_ann)), , drop = FALSE]

xp_bre_ann <- xp_ann[xp_ann$comparison == "breeding_improved_vs_landrace", , drop = FALSE]
xp_bre_ann <- xp_bre_ann[order(-xp_bre_ann$value), ][1:min(4, nrow(xp_bre_ann)), , drop = FALSE]

bal_cols <- bal_ann$color
bal_cols[is.na(bal_cols)] <- "#888888"

vz <- sort(unique(bal_ann$pos_mb))

# =========================
# Plot function
# =========================
plot_track <- function(df,
                       threshold,
                       ann_df = NULL,
                       v_guides = NULL,
                       y_min_data = NULL,
                       show_x = FALSE,
                       drop_zero_points = FALSE,
                       panel_name = "panel",
                       bottom_pad_frac = 0.06) {
  
  df_plot <- df
  
  if (drop_zero_points) {
    df_plot <- df_plot[is.finite(df_plot$value) & df_plot$value > 0, , drop = FALSE]
  }
  
  yvals <- df_plot$value
  if (!is.null(ann_df) && nrow(ann_df) > 0) {
    yvals <- c(yvals, ann_df$value)
  }
  
  ymax_data <- max(yvals, na.rm = TRUE)
  ymin_data <- min(yvals, na.rm = TRUE)
  
  if (!is.null(y_min_data)) {
    ymin_data <- y_min_data
  }
  
  if (ymax_data == ymin_data) ymax_data <- ymax_data + 1
  
  yrange <- ymax_data - ymin_data
  if (!is.finite(yrange) || yrange <= 0) yrange <- max(1, ymax_data)
  
  ymin_plot <- ymin_data - bottom_pad_frac * yrange
  
  if (!is.null(fixed_ticks[[panel_name]])) {
    yticks <- fixed_ticks[[panel_name]]
  } else {
    yticks <- seq(ymin_data, ymax_data, length.out = 4)
  }
  
  ymax_plot <- max(max(yticks), ymax_data) + 0.03 * yrange
  tick_store[[panel_name]] <<- yticks
  
  plot(df_plot$pos_mb, df_plot$value,
       col  = col_fn(df_plot$chr),
       pch  = 16,
       cex  = 0.12,
       bty  = "l",
       xlim = xlim_all,
       ylim = c(ymin_plot, ymax_plot),
       axes = FALSE,
       xlab = "",
       ylab = "",
       main = "",
       xaxs = "i",
       yaxs = "i")
  
  axis(2, at = yticks, las = 2, tck = -0.018, labels = FALSE)
  
  if (show_x) {
    axis(1,
         at = chr_centers_mb,
         labels = FALSE,
         tck = -0.03,
         lwd = 1,
         lwd.ticks = 1)
  }
  
  abline(h = threshold, col = "#FF0000", lty = 2, lwd = 0.8)
  
  if (!is.null(ann_df) && nrow(ann_df) > 0) {
    points(ann_df$pos_mb, ann_df$value,
           pch = 21,
           bg  = ann_df$color,
           col = NA,
           cex = 0.58,
           lwd = 0)
  }
  
  if (!is.null(v_guides) && length(v_guides) > 0) {
    abline(v = v_guides,
           col = "grey75",
           lty = 3,
           lwd = 0.5)
  }
}

# =========================
# Write TIFF
# =========================
dir.create(dirname(out_tiff), showWarnings = FALSE, recursive = TRUE)

tiff(out_tiff,
     width = 8,
     height = 7,
     units = "in",
     res = 600,
     compression = "lzw",
     type = "cairo",
     bg = "white")

par(mfrow = c(7, 1),
    mar = c(0.04, 1.15, 0.04, 0.12),
    oma = c(0.08, 0.10, 0.08, 0.08))

# 1. Wild B2
plot_track(
  wild_df,
  threshold = quantile(wild_df$value, 0.99, na.rm = TRUE),
  ann_df = data.frame(
    pos_mb = bal_ann$pos_mb,
    value  = bal_ann$wild,
    color  = bal_cols,
    stringsAsFactors = FALSE
  ),
  v_guides = vz,
  y_min_data = 0,
  show_x = FALSE,
  drop_zero_points = TRUE,
  panel_name = "wild_B2"
)

# 2. Landrace B2
plot_track(
  land_df,
  threshold = quantile(land_df$value, 0.99, na.rm = TRUE),
  ann_df = la_unique,
  v_guides = vz,
  y_min_data = 0,
  show_x = FALSE,
  drop_zero_points = TRUE,
  panel_name = "landrace_B2"
)

# 3. Improved B2
plot_track(
  imp_df,
  threshold = quantile(imp_df$value, 0.99, na.rm = TRUE),
  ann_df = imp_unique,
  v_guides = vz,
  y_min_data = 0,
  show_x = FALSE,
  drop_zero_points = TRUE,
  panel_name = "improved_B2"
)

# 4. FST wild vs landrace
plot_track(
  fst_wl_df,
  threshold = quantile(fst_wl_df$value, 0.99, na.rm = TRUE),
  ann_df = transform(fst_wl_ann, color = "#D62728"),
  show_x = FALSE,
  panel_name = "fst_wild_landrace"
)

# 5. XPCLR domestication
plot_track(
  xp_dom_df,
  threshold = quantile(xp_dom_df$value, 0.99, na.rm = TRUE),
  ann_df = transform(xp_dom_ann, color = "#D62728"),
  show_x = FALSE,
  panel_name = "xpclr_domestication"
)

# 6. FST landrace vs improved
plot_track(
  fst_li_df,
  threshold = quantile(fst_li_df$value, 0.99, na.rm = TRUE),
  ann_df = transform(fst_li_ann, color = "#D62728"),
  show_x = FALSE,
  panel_name = "fst_landrace_improved"
)

# 7. XPCLR breeding
plot_track(
  xp_bre_df,
  threshold = quantile(xp_bre_df$value, 0.99, na.rm = TRUE),
  ann_df = transform(xp_bre_ann, color = "#D62728"),
  show_x = TRUE,
  panel_name = "xpclr_breeding"
)

dev.off()

cat("Saved interactive final TIFF to:\n", out_tiff, "\n\n")

cat("Y-axis tick values used for each panel:\n")
for (nm in names(tick_store)) {
  cat("\n", nm, ": ", paste(format(round(tick_store[[nm]], 3), nsmall = 3), collapse = ", "), sep = "")
}
cat("\n")