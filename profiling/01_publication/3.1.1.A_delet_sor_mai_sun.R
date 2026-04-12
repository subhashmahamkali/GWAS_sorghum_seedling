#!/usr/bin/env Rscript
# =============================================================================
# Genome-Wide Absolute Genetic Load (Top 5%) - 3 Species Violin Plot
# Uses vioplot as suggested by advisor
# Sorghum:   wild (n=49) / landrace (n=101) / improved (n=128)
# Maize:     wild (n=20) / landrace (n=26)  / improved (n=46, no B73)
# Sunflower: wild (n=26) / landrace (n=18)  / improved (n=288)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(vioplot)
})

# --------------------------------------------------------------------------
# 0. Paths
# --------------------------------------------------------------------------
repo_root <- "/Users/subhashmahamkali/Documents/gwas_sap"
load_dir  <- file.path(repo_root,
                       "graphs/01_publication/3.genetic_load/genome_level_1_5_10_violin_custom_palette")
out_dir   <- load_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Read files
# --------------------------------------------------------------------------
sorg_all  <- fread(file.path(load_dir, "sorghum/sorghum_genome_level_by_sample_1_5_10.tsv"))
sun_all   <- fread(file.path(load_dir, "sunflower/sunflower_genome_level_by_sample_1_5_10.tsv"))
maize_all <- fread(file.path(load_dir, "maize/maize_5pct_wild_lan_imp_noB73.tsv"))

# --------------------------------------------------------------------------
# 2. Filter to 5% and extract per-group vectors
# --------------------------------------------------------------------------
sorg5 <- sorg_all[cutoff_lab == "5%"]
sun5  <- sun_all [cutoff_lab == "5%"]

sorg_wild  <- sorg5[group == "wild",     Abs_genetic_load]
sorg_land  <- sorg5[group == "landrace", Abs_genetic_load]
sorg_imp   <- sorg5[group == "improved", Abs_genetic_load]

sun_wild   <- sun5[group == "wild",     Abs_genetic_load]
sun_land   <- sun5[group == "landrace", Abs_genetic_load]
sun_imp    <- sun5[group == "improved", Abs_genetic_load]

maize_wild <- maize_all[group == "wild",     Abs_genetic_load]
maize_land <- maize_all[group == "landrace", Abs_genetic_load]
maize_imp  <- maize_all[group == "improved", Abs_genetic_load]

# --------------------------------------------------------------------------
# 3. Colours
# --------------------------------------------------------------------------
col_wild <- adjustcolor("#4DAF4A", alpha.f = 0.85)
col_land <- adjustcolor("#377EB8", alpha.f = 0.85)
col_imp  <- adjustcolor("#984EA3", alpha.f = 0.85)
cols     <- c(col_wild, col_land, col_imp)

# --------------------------------------------------------------------------
# 4. Y-axis label formatter: values >= 1M → "1.5M", else "300K" etc.
# --------------------------------------------------------------------------
#fmt_axis <- function(x) {
  #ifelse(abs(x) >= 1e6,
         #paste0(formatC(x / 1e6, format = "fg", digits = 2), "M"),
         #ifelse(abs(x) >= 1e3,
                #paste0(formatC(x / 1e3, format = "fg", digits = 3), "K"),
                #as.character(x)))
#}
fmt_axis <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

# --------------------------------------------------------------------------
# 5. Significance bracket helpers
# --------------------------------------------------------------------------
draw_bracket <- function(x1, x2, y, p_val, tick, label_offset) {
  sig <- if      (p_val < 0.0001) "****"
  else if (p_val < 0.001)  "***"
  else if (p_val < 0.01)   "**"
  else if (p_val < 0.05)   "*"
  else                     "ns"
  if (sig == "ns") return(invisible(NULL))
  segments(x1, y, x2, y, xpd = NA)
  segments(x1, y - tick, x1, y, xpd = NA)
  segments(x2, y - tick, x2, y, xpd = NA)
  text((x1 + x2) / 2, y + label_offset, labels = sig, cex = 1.0, xpd = NA)
}

add_wilcox_brackets <- function(g_list, x_pos, y_top, step) {
  pairs <- list(c(1,2), c(1,3), c(2,3))
  tick  <- step * 0.20
  loff  <- step * 0.22
  for (k in seq_along(pairs)) {
    i <- pairs[[k]][1]; j <- pairs[[k]][2]
    wt <- suppressWarnings(wilcox.test(g_list[[i]], g_list[[j]]))
    draw_bracket(x_pos[i], x_pos[j],
                 y            = y_top + (k - 1) * step,
                 p_val        = wt$p.value,
                 tick         = tick,
                 label_offset = loff)
  }
}

# --------------------------------------------------------------------------
# 6. Plotting function
# --------------------------------------------------------------------------
make_plot <- function() {
  par(mfrow = c(1, 3),
      mar   = c(3, 6, 3.5, 0.5),
      oma   = c(0, 0, 3, 0))
  
  # ── Panel A: Sorghum ────────────────────────────────────────────────────
  all_s  <- c(sorg_wild, sorg_land, sorg_imp)
  ylo_s  <- min(all_s) * 0.88
  step_s <- (max(all_s) - min(all_s)) * 0.11
  yhi_s  <- max(all_s) + step_s * 1.2
  
  vioplot(sorg_wild, sorg_land, sorg_imp,
          names   = c("wild", "landrace", "improved"),
          col     = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          ylim    = c(ylo_s, yhi_s),
          las = 1, tck = -0.03,
          ylab = "Absolute genetic load",
          main = "Sorghum",
          cex.lab = 1.05, cex.axis = 0.85, cex.main = 1.2, font.main = 2,
          yaxt = "n")
  grid(nx = NA, ny = NULL, col = "gray70", lty = "solid", lwd = 0.6)
  vioplot(sorg_wild, sorg_land, sorg_imp,
          col = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          add = TRUE)
  
  ticks_s <- pretty(c(ylo_s, max(all_s)), n = 5)
  axis(2, at = ticks_s, labels = fmt_axis(ticks_s), las = 1, cex.axis = 0.88)
  
  add_wilcox_brackets(list(sorg_wild, sorg_land, sorg_imp),
                      x_pos = 1:3, y_top = max(all_s) * 1.04, step = step_s)
  
  # ── Panel B: Maize ──────────────────────────────────────────────────────
  all_m  <- c(maize_wild, maize_land, maize_imp)
  ylo_m  <- min(all_m) * 0.75
  step_m <- (max(all_m) - min(all_m)) * 0.11
  yhi_m  <- max(all_m) + step_m * 1.2
  
  vioplot(maize_wild, maize_land, maize_imp,
          names   = c("wild", "landrace", "improved"),
          col     = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          ylim    = c(ylo_m, yhi_m),
          las = 1, tck = -0.03,
          ylab = "",
          main = "Maize",
          cex.lab = 1.05, cex.axis = 0.85, cex.main = 1.2, font.main = 2,
          yaxt = "n")
  grid(nx = NA, ny = NULL, col = "gray70", lty = "solid", lwd = 0.6)
  vioplot(maize_wild, maize_land, maize_imp,
          col = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          add = TRUE)
  
  ticks_m <- pretty(c(ylo_m, max(all_m)), n = 5)
  axis(2, at = ticks_m, labels = fmt_axis(ticks_m), las = 1, cex.axis = 0.88)
  
  add_wilcox_brackets(list(maize_wild, maize_land, maize_imp),
                      x_pos = 1:3, y_top = max(all_m) * 1.04, step = step_m)
  
  # ── Panel C: Sunflower ──────────────────────────────────────────────────
  all_sf  <- c(sun_wild, sun_land, sun_imp)
  ylo_sf  <- 0#min(all_sf) * 0.75
  step_sf <- (max(all_sf) - min(all_sf)) * 0.11
  yhi_sf  <- max(all_sf) + step_sf * 1.2
  
  vioplot(sun_wild, sun_land, sun_imp,
          names   = c("wild", "landrace", "improved"),
          col     = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          ylim    = c(ylo_sf, yhi_sf),
          las = 1, tck = -0.03,
          ylab = "",
          main = "Sunflower",
          cex.lab = 1.05, cex.axis = 0.85, cex.main = 1.2, font.main = 2,
          yaxt = "n")
  grid(nx = NA, ny = NULL, col = "gray70", lty = "solid", lwd = 0.6)
  vioplot(sun_wild, sun_land, sun_imp,
          col = cols, border = cols,
          rectCol = "white", lineCol = "black", colMed = "black",
          add = TRUE)
  
  ticks_sf <- pretty(c(0, max(all_sf)), n = 5)
  axis(2, at = ticks_sf, labels = fmt_axis(ticks_sf), las = 1, cex.axis = 0.88)
  
  add_wilcox_brackets(list(sun_wild, sun_land, sun_imp),
                      x_pos = 1:3, y_top = max(all_sf) * 1.04, step = step_sf)
  # ── Overall title ────────────────────────────────────────────────────────
  mtext("Genome-Wide Absolute Genetic Load (Top 5%)",
        outer = TRUE, side = 3, line = 1, cex = 1.3, font = 2)
}

# --------------------------------------------------------------------------
# 7. Save PDF + PNG
# --------------------------------------------------------------------------
pdf(file.path(out_dir, "s_m_s_5pct_3panel_vioplot.pdf"),
    width = 13, height = 5.5)
make_plot()
dev.off()

png(file.path(out_dir, "s_m_s_5pct_3panel_vioplot.png"),
    width = 13, height = 5.5, units = "in", res = 350, bg = "white")
make_plot()
dev.off()

message("Done. Outputs saved to: ", out_dir)