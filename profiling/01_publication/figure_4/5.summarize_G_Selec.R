#!/usr/bin/env Rscript
# ==============================================================================
# GWAS loci and selection overlap by curated trait category
# Updated to follow actual GWAS trait parsing logic
#
# KEY FEATURES
# 1. Uses pattern-based curated trait mapping from Trait_full
# 2. Derives nitrogen condition from Trait_full (HN / LN / NR)
# 3. Uses only selected curated traits
# 4. Builds GWAS loci as +/- 40 kb around significant SNPs, then merges overlaps
# 5. Counts overlap separately for:
#      Balancing: wild, landrace, improved
#      Positive : FST wild-landrace, FST landrace-improved,
#                 XPCLR domestication, XPCLR breeding
# 6. Summarizes both:
#      - sum = overlap events
#      - any = unique loci overlapping at least one scan
# 7. Produces MBE-style plot with legend
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

OUTDIR <- file.path(repo_root, "largedata", "figure_4", "panel_A_curated_traits_MBE_v2")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PARAMETERS
# ==============================================================================
GWAS_THRESHOLD <- 5
TOP_PERCENTILE <- 0.01

GWAS_FLANK    <- 40000   # +/- 40 kb
B2_MERGE_DIST <- 10000
B2_BUFFER     <- 5000

# ==============================================================================
# HELPERS
# ==============================================================================
merge_intervals <- function(dt) {
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
    ends <- integer()
    
    if (.N > 1) {
      for (i in 2:.N) {
        if (start[i] <= cur_e + 1L) {
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
    ends <- integer()
    
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

# expected columns from your original summary script
setnames(gwas, c("CHROM", "POS", "Pvalue", "Trait", "Condition", "Trait_type"))

gwas[, CHROM := as.integer(CHROM)]
gwas[, POS   := as.integer(POS)]
gwas[, logp  := -log10(Pvalue)]

# keep original trait string
gwas[, Trait_full := Trait]

# ==============================================================================
# CURATED TRAIT MAPPING
# Based on your actual GWAS parsing code
# ==============================================================================
trait_map <- data.table(
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
  category = c(
    "Architecture",
    "Architecture",
    "Architecture",
    "Architecture",
    rep("Developmental", 10),
    rep("Panicle", 6)
  )
)

gwas[, Trait_clean := NA_character_]
gwas[, category := NA_character_]

for (i in seq_len(nrow(trait_map))) {
  hit <- grepl(trait_map$pattern[i], gwas$Trait_full, fixed = TRUE)
  gwas$Trait_clean[hit] <- trait_map$Trait_clean[i]
  gwas$category[hit]    <- trait_map$category[i]
}

# keep only curated traits
gwas2 <- gwas[!is.na(Trait_clean)]

# derive nitrogen treatment from Trait_full, matching your actual code logic
gwas2[, Nitrogen_Treatment := NA_character_]
gwas2[grepl("\\.HN\\.|^HN\\.", Trait_full), Nitrogen_Treatment := "HN"]
gwas2[grepl("\\.LN\\.|^LN\\.", Trait_full), Nitrogen_Treatment := "LN"]
gwas2[grepl("\\.NR$|\\.NR\\.|^NR\\.|NR$", Trait_full), Nitrogen_Treatment := "NR"]

gwas2 <- gwas2[!is.na(Nitrogen_Treatment)]

if (nrow(gwas2) == 0) {
  stop("No GWAS rows remained after curated trait mapping and nitrogen parsing.")
}

cat("\nMatched curated traits:\n")
print(unique(gwas2[, .(Trait_full, Trait_clean, category)])[order(category, Trait_clean)])

# significant SNPs only
gwas_sig <- gwas2[
  is.finite(logp) &
    logp > GWAS_THRESHOLD
]

if (nrow(gwas_sig) == 0) {
  stop("No significant GWAS SNPs remained after filtering.")
}

cat("\nSignificant GWAS SNPs by category and nitrogen condition:\n")
print(gwas_sig[, .N, by = .(category, Nitrogen_Treatment)][order(category, Nitrogen_Treatment)])

fwrite(gwas_sig, file.path(OUTDIR, "gwas_significant_curated_traits.tsv"), sep = "\t")

# ==============================================================================
# BUILD GWAS LOCI
# ==============================================================================
cat("\nBuilding GWAS loci using +/-", GWAS_FLANK/1000, "kb windows...\n")

loci_list <- list()
idx <- 0L

for (cond in c("HN", "LN", "NR")) {
  for (catg in c("Architecture", "Developmental", "Panicle")) {
    
    sub <- gwas_sig[Nitrogen_Treatment == cond & category == catg]
    if (nrow(sub) == 0) next
    
    loci <- make_gwas_loci_pm40kb(sub$CHROM, sub$POS, flank = GWAS_FLANK)
    if (nrow(loci) == 0) next
    
    loci[, `:=`(Condition = cond, category = catg)]
    setcolorder(loci, c("Condition", "category", "locus_id", "chr", "start", "end"))
    
    idx <- idx + 1L
    loci_list[[idx]] <- loci
    
    cat("  ", cond, "x", catg, ":",
        nrow(unique(sub[, .(CHROM, POS)])), "significant SNPs ->",
        nrow(loci), "merged loci\n")
  }
}

if (length(loci_list) == 0) {
  stop("No GWAS loci were created. Check trait parsing and nitrogen parsing.")
}

all_loci <- rbindlist(loci_list, use.names = TRUE, fill = TRUE)

cat("\nGWAS loci by category and condition:\n")
print(all_loci[, .N, by = .(category, Condition)][order(category, Condition)])

fwrite(all_loci, file.path(OUTDIR, "gwas_loci_pm40kb_merged_curated_traits.tsv"), sep = "\t")

# ==============================================================================
# LOAD BALANCING PEAKS
# ==============================================================================
cat("\nLoading balancing peaks...\n")

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
  merge_snps_to_intervals(out_top, merge_dist = B2_MERGE_DIST, buffer = B2_BUFFER)
}

b2_wild <- load_b2(b2_wild_path, "wild")
b2_land <- load_b2(b2_land_path, "landrace")
b2_imp  <- load_b2(b2_imp_path, "improved")

cat("  B2 wild peaks     :", nrow(b2_wild), "\n")
cat("  B2 landrace peaks :", nrow(b2_land), "\n")
cat("  B2 improved peaks :", nrow(b2_imp), "\n")

# ==============================================================================
# LOAD POSITIVE SELECTION PEAKS
# ==============================================================================
cat("\nLoading positive-selection peaks...\n")

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

cat("  FST wild-landrace peaks     :", nrow(fst_wl), "\n")
cat("  FST landrace-improved peaks :", nrow(fst_li), "\n")
cat("  XPCLR domestication peaks   :", nrow(xp_dom), "\n")
cat("  XPCLR breeding peaks        :", nrow(xp_bre), "\n")

# ==============================================================================
# OVERLAPS
# ==============================================================================
cat("\nComputing locus overlaps...\n")

all_loci[, bal_wild := interval_overlap_any(.SD, b2_wild), .SDcols = c("chr", "start", "end")]
all_loci[, bal_land := interval_overlap_any(.SD, b2_land), .SDcols = c("chr", "start", "end")]
all_loci[, bal_imp  := interval_overlap_any(.SD, b2_imp),  .SDcols = c("chr", "start", "end")]

all_loci[, pos_fst_wl := interval_overlap_any(.SD, fst_wl), .SDcols = c("chr", "start", "end")]
all_loci[, pos_fst_li := interval_overlap_any(.SD, fst_li), .SDcols = c("chr", "start", "end")]
all_loci[, pos_xp_dom := interval_overlap_any(.SD, xp_dom), .SDcols = c("chr", "start", "end")]
all_loci[, pos_xp_bre := interval_overlap_any(.SD, xp_bre), .SDcols = c("chr", "start", "end")]

all_loci[, bal_any := bal_wild | bal_land | bal_imp]
all_loci[, pos_any := pos_fst_wl | pos_fst_li | pos_xp_dom | pos_xp_bre]

all_loci[, bal_sum := as.integer(bal_wild) + as.integer(bal_land) + as.integer(bal_imp)]
all_loci[, pos_sum := as.integer(pos_fst_wl) + as.integer(pos_fst_li) +
           as.integer(pos_xp_dom) + as.integer(pos_xp_bre)]

fwrite(all_loci, file.path(OUTDIR, "gwas_loci_with_overlap_flags_curated_traits.tsv"), sep = "\t")

# ==============================================================================
# SUMMARY
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

cat("\nSummary table:\n")
print(results)

results_pct <- copy(results)
for (cc in c("n_bal_wild","n_bal_land","n_bal_imp","n_bal_sum","n_bal_any",
             "n_pos_fst_wl","n_pos_fst_li","n_pos_xp_dom","n_pos_xp_bre","n_pos_sum","n_pos_any")) {
  results_pct[, paste0("pct_", cc) := round(100 * get(cc) / n_loci, 1)]
}

fwrite(results,     file.path(OUTDIR, "summary_counts_curated_traits.tsv"), sep = "\t")
fwrite(results_pct, file.path(OUTDIR, "summary_percentages_curated_traits.tsv"), sep = "\t")

# ==============================================================================
# MAIN MBE-STYLE PLOT
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
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
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
      override.aes = list(linewidth = 1.3, shape = 16, size = 2.6),
      nrow = 1
    )
  )

ggsave(
  file.path(OUTDIR, "plot_unique_any_bal_pos_MBE_curated_traits.pdf"),
  p_mbe,
  width = 10.5,
  height = 5.8,
  dpi = 300,
  useDingbats = FALSE
)

ggsave(
  file.path(OUTDIR, "plot_unique_any_bal_pos_MBE_curated_traits.png"),
  p_mbe,
  width = 10.5,
  height = 5.8,
  dpi = 600
)

cat("\nSaved outputs in:\n", OUTDIR, "\n")
cat("\nDone.\n")