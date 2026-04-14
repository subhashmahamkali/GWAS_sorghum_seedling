#!/usr/bin/env Rscript
# ==============================================================================
# GWAS loci and selection overlap by trait category
# Publication-style summary figure (MBE style)
#
# LOGIC
# 1. Load significant GWAS SNPs
# 2. Expand each SNP to +/- 40 kb
# 3. Merge overlapping SNP windows into non-redundant GWAS loci
# 4. Build selection peaks separately:
#      Balancing: wild, landrace, improved
#      Positive : FST wild-landrace, FST landrace-improved,
#                 XPCLR domestication, XPCLR breeding
# 5. Count overlap of GWAS loci with each scan separately
# 6. Summarize:
#      - separate counts
#      - unique balancing overlap (any balancing scan)
#      - unique positive overlap (any positive scan)
# 7. Make MBE-style plot
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ==============================================================================
# PATHS
# ==============================================================================
repo_root <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling"

gwas_path    <- file.path(repo_root, "largedata", "All_GWAS_sorg.txt")
b2_wild_path <- file.path(repo_root, "largedata", "wild_merged_chr.B2.txt")
b2_land_path <- file.path(repo_root, "largedata", "Sorghum_landrace.AGPV5.B2_stat.txt")
b2_imp_path  <- file.path(repo_root, "largedata", "imp_merged_chr.B2.txt")

fst_wl_path  <- file.path(repo_root, "largedata", "figure_2", "sorg_wild_vs_landrace.windowed.weir.fst")
fst_li_path  <- file.path(repo_root, "largedata", "figure_2", "landrace_vs_improved.windowed.weir.fst")
xp_dom_path  <- file.path(repo_root, "largedata", "figure_2", "domestication_landrace_vs_wild.all_chr.tsv")
xp_bre_path  <- file.path(repo_root, "largedata", "figure_2", "breeding_improved_vs_landrace.all_chr.tsv")

OUTDIR <- file.path(repo_root, "largedata", "figure_4", "panel_A_summary_split_selection_MBE")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PARAMETERS
# ==============================================================================
GWAS_THRESHOLD <- 5
TOP_PERCENTILE <- 0.01

GWAS_FLANK    <- 40000   # +/- 40 kb around each significant GWAS SNP
B2_MERGE_DIST <- 10000   # merge nearby B2 top SNPs within each population
B2_BUFFER     <- 5000    # extend merged B2 regions by 5 kb

# ==============================================================================
# HELPERS
# ==============================================================================

merge_intervals <- function(dt) {
  # dt must contain: chr, start, end
  if (nrow(dt) == 0) {
    return(data.table(chr = integer(), start = integer(), end = integer()))
  }
  
  dt <- copy(dt)
  dt[, chr := as.integer(chr)]
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  setorder(dt, chr, start, end)
  
  out <- dt[, {
    cur_s <- start[1]
    cur_e <- end[1]
    starts <- integer()
    ends   <- integer()
    
    if (.N > 1) {
      for (i in 2:.N) {
        if (start[i] <= (cur_e + 1L)) {
          cur_e <- max(cur_e, end[i])
        } else {
          starts <- c(starts, cur_s)
          ends   <- c(ends, cur_e)
          cur_s  <- start[i]
          cur_e  <- end[i]
        }
      }
    }
    
    starts <- c(starts, cur_s)
    ends   <- c(ends, cur_e)
    .(start = starts, end = ends)
  }, by = chr]
  
  out[]
}

merge_snps_to_intervals <- function(dt, merge_dist = 10000, buffer = 5000) {
  # dt must contain: chr, pos
  if (nrow(dt) == 0) {
    return(data.table(chr = integer(), start = integer(), end = integer()))
  }
  
  dt <- copy(dt)
  dt[, chr := as.integer(chr)]
  dt[, pos := as.integer(pos)]
  setorder(dt, chr, pos)
  
  out <- dt[, {
    s <- pos[1]
    e <- pos[1]
    starts <- integer()
    ends   <- integer()
    
    if (.N > 1) {
      for (i in 2:.N) {
        if ((pos[i] - e) <= merge_dist) {
          e <- pos[i]
        } else {
          starts <- c(starts, s)
          ends   <- c(ends, e)
          s <- pos[i]
          e <- pos[i]
        }
      }
    }
    
    starts <- c(starts, s)
    ends   <- c(ends, e)
    .(
      start = pmax(1L, starts - buffer),
      end   = ends + buffer
    )
  }, by = chr]
  
  out[]
}

make_gwas_loci_pm40kb <- function(chr_vec, pos_vec, flank = 40000) {
  dt <- data.table(
    chr = as.integer(chr_vec),
    pos = as.integer(pos_vec)
  )
  dt <- unique(dt)
  
  dt[, start := pmax(1L, pos - flank)]
  dt[, end   := pos + flank]
  
  loci <- merge_intervals(dt[, .(chr, start, end)])
  loci[, locus_id := paste0("L", seq_len(.N))]
  setcolorder(loci, c("locus_id", "chr", "start", "end"))
  loci[]
}

interval_overlap_any <- function(query_dt, peaks_dt) {
  # query_dt: chr, start, end
  # peaks_dt: chr, start, end
  if (nrow(query_dt) == 0) return(logical())
  if (nrow(peaks_dt) == 0) return(rep(FALSE, nrow(query_dt)))
  
  q <- copy(query_dt)
  p <- copy(peaks_dt)
  
  q[, qid := .I]
  
  setkey(q, chr, start, end)
  setkey(p, chr, start, end)
  
  ov <- foverlaps(q, p, type = "any", nomatch = 0L)
  hit_ids <- unique(ov$qid)
  
  out <- rep(FALSE, nrow(q))
  out[hit_ids] <- TRUE
  out
}

# ==============================================================================
# LOAD GWAS
# ==============================================================================
cat("Loading GWAS data...\n")

gwas <- fread(gwas_path)
setnames(gwas, c("CHROM", "POS", "Pvalue", "Trait", "Condition", "Trait_type"))

gwas[, CHROM := as.integer(CHROM)]
gwas[, POS   := as.integer(POS)]
gwas[, logp  := -log10(Pvalue)]

# Trait category mapping
gwas[, category := fcase(
  grepl("plantHeight", Trait, ignore.case = TRUE) & !grepl("panicle", Trait, ignore.case = TRUE), "Architecture",
  grepl("tillersPerPlant", Trait, ignore.case = TRUE), "Architecture",
  grepl("branchInternodeLength", Trait, ignore.case = TRUE), "Architecture",
  grepl("primaryBranchNumber", Trait, ignore.case = TRUE), "Architecture",
  
  grepl("daysToFlower", Trait, ignore.case = TRUE), "Developmental",
  grepl("extantLeafNumber", Trait, ignore.case = TRUE), "Developmental",
  grepl("flagLeaf", Trait, ignore.case = TRUE), "Developmental",
  grepl("thirdLeaf", Trait, ignore.case = TRUE), "Developmental",
  grepl("leafAngle", Trait, ignore.case = TRUE), "Developmental",
  grepl("medianLeafAngle", Trait, ignore.case = TRUE), "Developmental",
  grepl("stemDiameter", Trait, ignore.case = TRUE), "Developmental",
  
  grepl("panicle", Trait, ignore.case = TRUE), "Panicle",
  grepl("rachis", Trait, ignore.case = TRUE), "Panicle",
  grepl("estimatedPlotYield", Trait, ignore.case = TRUE), "Panicle",
  grepl("panicleGrainWeight", Trait, ignore.case = TRUE), "Panicle",
  grepl("seedMass", Trait, ignore.case = TRUE), "Panicle",
  grepl("percent", Trait, ignore.case = TRUE), "Panicle",
  
  default = "Other"
)]

gwas_sig <- gwas[
  is.finite(logp) &
    logp > GWAS_THRESHOLD &
    category != "Other"
]

cat("\nSignificant GWAS SNPs by condition and category:\n")
print(gwas_sig[, .N, by = .(Condition, category)][order(category, Condition)])

# ==============================================================================
# BUILD GWAS LOCI: +/- 40 kb, merge overlaps
# ==============================================================================
cat("\nBuilding GWAS loci using +/-", GWAS_FLANK / 1000, "kb and merging overlaps...\n")

loci_list <- list()
idx <- 0L

for (cond in c("HN", "LN", "NR")) {
  for (catg in c("Architecture", "Developmental", "Panicle")) {
    
    sub <- gwas_sig[Condition == cond & category == catg]
    if (nrow(sub) == 0) next
    
    loci <- make_gwas_loci_pm40kb(sub$CHROM, sub$POS, flank = GWAS_FLANK)
    loci[, `:=`(Condition = cond, category = catg)]
    setcolorder(loci, c("Condition", "category", "locus_id", "chr", "start", "end"))
    
    idx <- idx + 1L
    loci_list[[idx]] <- loci
    
    cat("  ", cond, "x", catg, ":",
        nrow(unique(sub[, .(CHROM, POS)])), "significant SNPs ->",
        nrow(loci), "merged loci\n")
  }
}

all_loci <- rbindlist(loci_list, use.names = TRUE, fill = TRUE)

cat("\nGWAS loci by condition and category:\n")
print(all_loci[, .N, by = .(Condition, category)][order(category, Condition)])

fwrite(all_loci, file.path(OUTDIR, "gwas_loci_pm40kb_merged.tsv"), sep = "\t")

# ==============================================================================
# LOAD BALANCING SELECTION PEAKS SEPARATELY
# ==============================================================================
cat("\nLoading balancing selection peaks...\n")

load_b2 <- function(path, pop) {
  d <- fread(path)
  
  if (pop == "landrace") {
    out <- data.table(
      chr   = as.integer(d[[1]]),
      pos   = as.integer(d[[2]]),
      score = as.numeric(d[[4]])
    )
  } else {
    out <- data.table(
      chr   = as.integer(d[[7]]),
      pos   = as.integer(d[[1]]),
      score = as.numeric(d[[3]])
    )
  }
  
  out <- out[is.finite(score) & score > 0]
  thr <- quantile(out$score, 1 - TOP_PERCENTILE, na.rm = TRUE)
  
  out_top <- out[score >= thr, .(chr, pos)]
  peaks <- merge_snps_to_intervals(out_top, merge_dist = B2_MERGE_DIST, buffer = B2_BUFFER)
  peaks[]
}

b2_wild <- load_b2(b2_wild_path, "wild")
b2_land <- load_b2(b2_land_path, "landrace")
b2_imp  <- load_b2(b2_imp_path, "improved")

cat("  B2 wild peaks     :", nrow(b2_wild), "\n")
cat("  B2 landrace peaks :", nrow(b2_land), "\n")
cat("  B2 improved peaks :", nrow(b2_imp), "\n")

# ==============================================================================
# LOAD POSITIVE SELECTION PEAKS SEPARATELY
# ==============================================================================
cat("\nLoading positive selection peaks...\n")

load_fst <- function(path) {
  d <- fread(path)
  d[, CHROM := as.integer(CHROM)]
  d <- d[is.finite(WEIGHTED_FST)]
  
  thr <- quantile(d$WEIGHTED_FST, 1 - TOP_PERCENTILE, na.rm = TRUE)
  
  peaks <- d[WEIGHTED_FST >= thr,
             .(chr = CHROM,
               start = as.integer(BIN_START),
               end   = as.integer(BIN_END))]
  
  merge_intervals(peaks)
}

load_xp <- function(path) {
  d <- fread(path)
  d[, chrom := as.integer(chrom)]
  d <- d[is.finite(xpclr) & xpclr > 0]
  
  thr <- quantile(d$xpclr, 1 - TOP_PERCENTILE, na.rm = TRUE)
  
  peaks <- d[xpclr >= thr,
             .(chr = chrom,
               start = as.integer(start),
               end   = as.integer(stop))]
  
  merge_intervals(peaks)
}

fst_wl <- load_fst(fst_wl_path)
fst_li <- load_fst(fst_li_path)
xp_dom <- load_xp(xp_dom_path)
xp_bre <- load_xp(xp_bre_path)

cat("  FST wild-landrace peaks       :", nrow(fst_wl), "\n")
cat("  FST landrace-improved peaks   :", nrow(fst_li), "\n")
cat("  XPCLR domestication peaks     :", nrow(xp_dom), "\n")
cat("  XPCLR breeding peaks          :", nrow(xp_bre), "\n")

# save peaks
fwrite(b2_wild, file.path(OUTDIR, "balancing_wild_peaks.tsv"), sep = "\t")
fwrite(b2_land, file.path(OUTDIR, "balancing_landrace_peaks.tsv"), sep = "\t")
fwrite(b2_imp,  file.path(OUTDIR, "balancing_improved_peaks.tsv"), sep = "\t")
fwrite(fst_wl,  file.path(OUTDIR, "positive_fst_wild_landrace_peaks.tsv"), sep = "\t")
fwrite(fst_li,  file.path(OUTDIR, "positive_fst_landrace_improved_peaks.tsv"), sep = "\t")
fwrite(xp_dom,  file.path(OUTDIR, "positive_xpclr_domestication_peaks.tsv"), sep = "\t")
fwrite(xp_bre,  file.path(OUTDIR, "positive_xpclr_breeding_peaks.tsv"), sep = "\t")

# ==============================================================================
# OVERLAP: separate scans first
# ==============================================================================
cat("\nComputing GWAS locus overlap with each scan separately...\n")

all_loci[, bal_wild := interval_overlap_any(.SD, b2_wild), .SDcols = c("chr", "start", "end")]
all_loci[, bal_land := interval_overlap_any(.SD, b2_land), .SDcols = c("chr", "start", "end")]
all_loci[, bal_imp  := interval_overlap_any(.SD, b2_imp),  .SDcols = c("chr", "start", "end")]

all_loci[, pos_fst_wl := interval_overlap_any(.SD, fst_wl), .SDcols = c("chr", "start", "end")]
all_loci[, pos_fst_li := interval_overlap_any(.SD, fst_li), .SDcols = c("chr", "start", "end")]
all_loci[, pos_xp_dom := interval_overlap_any(.SD, xp_dom), .SDcols = c("chr", "start", "end")]
all_loci[, pos_xp_bre := interval_overlap_any(.SD, xp_bre), .SDcols = c("chr", "start", "end")]

# unique counts
all_loci[, bal_any := bal_wild | bal_land | bal_imp]
all_loci[, pos_any := pos_fst_wl | pos_fst_li | pos_xp_dom | pos_xp_bre]

# summed events
all_loci[, bal_sum := as.integer(bal_wild) + as.integer(bal_land) + as.integer(bal_imp)]
all_loci[, pos_sum := as.integer(pos_fst_wl) + as.integer(pos_fst_li) +
           as.integer(pos_xp_dom) + as.integer(pos_xp_bre)]

fwrite(all_loci, file.path(OUTDIR, "gwas_loci_with_separate_overlap_flags.tsv"), sep = "\t")

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================
results <- all_loci[, .(
  n_loci = .N,
  
  n_bal_wild = sum(bal_wild),
  n_bal_land = sum(bal_land),
  n_bal_imp  = sum(bal_imp),
  n_bal_sum  = sum(bal_sum),
  n_bal_any  = sum(bal_any),
  
  n_pos_fst_wl = sum(pos_fst_wl),
  n_pos_fst_li = sum(pos_fst_li),
  n_pos_xp_dom = sum(pos_xp_dom),
  n_pos_xp_bre = sum(pos_xp_bre),
  n_pos_sum    = sum(pos_sum),
  n_pos_any    = sum(pos_any)
), by = .(Condition, category)]

results[, Condition := factor(Condition, levels = c("HN", "LN", "NR"))]
results[, category  := factor(category, levels = c("Architecture", "Developmental", "Panicle"))]
setorder(results, category, Condition)

cat("\nSummary results:\n")
print(results)

fwrite(results, file.path(OUTDIR, "independent_loci_separate_and_combined_overlap.tsv"), sep = "\t")

# ==============================================================================
# LONG TABLES
# ==============================================================================
sep_long <- melt(
  results,
  id.vars = c("Condition", "category", "n_loci"),
  measure.vars = c(
    "n_bal_wild", "n_bal_land", "n_bal_imp",
    "n_pos_fst_wl", "n_pos_fst_li", "n_pos_xp_dom", "n_pos_xp_bre"
  ),
  variable.name = "scan",
  value.name = "count"
)

sep_long[, selection_class := fifelse(grepl("^n_bal_", scan), "Balancing", "Positive")]
sep_long[, scan_label := fcase(
  scan == "n_bal_wild",   "B2 wild",
  scan == "n_bal_land",   "B2 landrace",
  scan == "n_bal_imp",    "B2 improved",
  scan == "n_pos_fst_wl", "FST wild-landrace",
  scan == "n_pos_fst_li", "FST landrace-improved",
  scan == "n_pos_xp_dom", "XPCLR domestication",
  scan == "n_pos_xp_bre", "XPCLR breeding"
)]

fwrite(sep_long, file.path(OUTDIR, "separate_overlap_counts_long.tsv"), sep = "\t")

results_pct <- copy(results)
for (cc in c("n_bal_wild","n_bal_land","n_bal_imp","n_bal_sum","n_bal_any",
             "n_pos_fst_wl","n_pos_fst_li","n_pos_xp_dom","n_pos_xp_bre","n_pos_sum","n_pos_any")) {
  results_pct[, paste0("pct_", cc) := round(100 * get(cc) / n_loci, 1)]
}
fwrite(results_pct, file.path(OUTDIR, "independent_loci_overlap_percentages.tsv"), sep = "\t")

# ==============================================================================
# MBE-STYLE PLOT
# Grey bars = total GWAS loci
# Orange    = unique balancing overlap (n_bal_any)
# Red       = unique positive overlap (n_pos_any)
# ==============================================================================

marker_dt <- melt(
  results,
  id.vars = c("Condition", "category", "n_loci"),
  measure.vars = c("n_bal_any", "n_pos_any"),
  variable.name = "selection",
  value.name = "count"
)

marker_dt[, selection := fcase(
  selection == "n_bal_any", "Balancing selection",
  selection == "n_pos_any", "Positive selection"
)]
marker_dt[, selection := factor(selection,
                                levels = c("Balancing selection", "Positive selection"))]

marker_dt[, Condition := factor(Condition, levels = c("HN", "LN", "NR"))]
marker_dt[, category  := factor(category, levels = c("Architecture", "Developmental", "Panicle"))]

ymax <- max(results$n_loci)
ybreaks <- seq(0, ceiling((ymax + 20) / 50) * 50, by = 50)

p_mbe <- ggplot() +
  geom_col(
    data = results,
    aes(x = Condition, y = n_loci),
    fill = "#B7C1C8",
    color = "#4D4D4D",
    width = 0.56,
    linewidth = 0.35
  ) +
  geom_errorbar(
    data = marker_dt,
    aes(x = Condition, ymin = count, ymax = count, color = selection),
    width = 0.36,
    linewidth = 1.15,
    show.legend = TRUE
  ) +
  geom_point(
    data = marker_dt,
    aes(x = Condition, y = count, color = selection),
    size = 2.2,
    stroke = 0,
    show.legend = TRUE
  ) +
  geom_text(
    data = results,
    aes(x = Condition, y = n_loci, label = n_loci),
    vjust = -0.38,
    size = 4.1,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = marker_dt,
    aes(x = Condition, y = count, label = count, color = selection),
    vjust = 1.75,
    size = 3.15,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~category, nrow = 1) +
  scale_color_manual(
    values = c(
      "Balancing selection" = "#E69F00",
      "Positive selection"  = "#CC3311"
    )
  ) +
  scale_y_continuous(
    breaks = ybreaks,
    limits = c(0, max(ybreaks)),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Nitrogen condition",
    y = "Number of GWAS loci",
    color = NULL,
    title = "GWAS loci overlapping balancing and positive selection"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 11.5, color = "black"),
    axis.text.y  = element_text(size = 11.5, color = "black"),
    axis.line    = element_line(color = "black", linewidth = 0.7),
    axis.ticks   = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.16, "cm"),
    
    strip.background = element_rect(fill = "#F2F2F2", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 12.5, color = "black"),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 11.5),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.45, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.3,
        shape = 16,
        size = 2.6
      ),
      nrow = 1
    )
  )

ggsave(
  file.path(OUTDIR, "plot_unique_any_bal_pos_MBE.pdf"),
  p_mbe,
  width = 10.5,
  height = 5.8,
  dpi = 300,
  useDingbats = FALSE
)

ggsave(
  file.path(OUTDIR, "plot_unique_any_bal_pos_MBE.png"),
  p_mbe,
  width = 10.5,
  height = 5.8,
  dpi = 600
)

# ==============================================================================
# OPTIONAL PLOT: separate scan counts
# ==============================================================================
sep_long[, scan_label := factor(
  scan_label,
  levels = c(
    "B2 wild", "B2 landrace", "B2 improved",
    "FST wild-landrace", "FST landrace-improved",
    "XPCLR domestication", "XPCLR breeding"
  )
)]

p_sep <- ggplot(sep_long, aes(x = Condition, y = count, fill = scan_label)) +
  geom_col(
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "grey25",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.82),
    vjust = -0.28,
    size = 2.7
  ) +
  facet_wrap(~category, nrow = 1) +
  scale_fill_manual(values = c(
    "B2 wild"               = "#F6C85F",
    "B2 landrace"           = "#E69F00",
    "B2 improved"           = "#C17C00",
    "FST wild-landrace"     = "#F28E8E",
    "FST landrace-improved" = "#E15759",
    "XPCLR domestication"   = "#C73E1D",
    "XPCLR breeding"        = "#8C1D18"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = "Nitrogen condition",
    y = "Number of overlapping GWAS loci",
    fill = NULL,
    title = "Separate overlap counts across balancing and positive selection scans"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "#F2F2F2", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank()
  )

ggsave(
  file.path(OUTDIR, "plot_separate_scans_grouped.pdf"),
  p_sep,
  width = 13.2,
  height = 5.6,
  dpi = 300,
  useDingbats = FALSE
)

# ==============================================================================
# DONE
# ==============================================================================
cat("\nSaved outputs in:\n", OUTDIR, "\n")
cat("\nMain MBE figure files:\n")
cat("  - plot_unique_any_bal_pos_MBE.pdf\n")
cat("  - plot_unique_any_bal_pos_MBE.png\n")
cat("\nDone.\n")