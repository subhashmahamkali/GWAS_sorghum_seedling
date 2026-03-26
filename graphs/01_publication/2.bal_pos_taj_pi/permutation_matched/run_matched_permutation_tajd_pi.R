#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(20260316)

get_env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (x == "") return(default)
  as.integer(x)
}

N_PERM <- get_env_int("N_PERM", 10000L)

BASE_DIR <- "/Users/subhashmahamkali/Documents/gwas_sap"
IN_DIR <- file.path(BASE_DIR, "data/2.1_TajimasD_PI")
OUT_DIR <- file.path(BASE_DIR, "graphs/01_publication/2.bal_pos_taj_pi/permutation_matched")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_tajd <- function(path, pop_label) {
  dt <- fread(
    path,
    col.names = c("chr", "start", "end", "stat", "set", "pop"),
    showProgress = FALSE
  )
  dt[, `:=`(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    stat = as.numeric(stat),
    set = as.character(set),
    pop = pop_label,
    metric = "TajimaD"
  )]
  dt[is.finite(stat) & set %in% c("PEAK", "BACKGROUND")]
}

read_pi <- function(path, pop_label) {
  dt <- fread(
    path,
    col.names = c("chr", "start", "end", "stat", "n_vars", "set", "pop"),
    showProgress = FALSE
  )
  dt[, `:=`(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    stat = as.numeric(stat),
    set = as.character(set),
    pop = pop_label,
    metric = "Pi"
  )]
  dt[is.finite(stat) & set %in% c("PEAK", "BACKGROUND")]
}

make_peak_regions <- function(dt) {
  pk <- copy(dt[set == "PEAK", .(chr, start, end, stat)])
  setorder(pk, chr, start, end)
  if (nrow(pk) == 0L) {
    return(data.table(
      chr = character(),
      grp = integer(),
      region_start = integer(),
      region_end = integer(),
      n_windows = integer(),
      region_stat = numeric(),
      region_id = character()
    ))
  }
  pk[, grp := cumsum(c(TRUE, start[-1] > end[-.N])), by = chr]
  reg <- pk[, .(
    region_start = min(start),
    region_end = max(end),
    n_windows = .N,
    region_stat = median(stat, na.rm = TRUE)
  ), by = .(chr, grp)]
  reg[, region_id := sprintf("%s_%04d", chr, seq_len(.N)), by = chr]
  reg[]
}

make_bg_struct <- function(dt) {
  bg <- copy(dt[set == "BACKGROUND", .(chr, start, end, stat)])
  setorder(bg, chr, start, end)
  if (nrow(bg) == 0L) return(list())
  bg[, run_id := cumsum(c(TRUE, start[-1] > end[-.N])), by = chr]
  by_chr <- split(bg, by = "chr", keep.by = FALSE)
  out <- vector("list", length(by_chr))
  names(out) <- names(by_chr)
  for (chr_name in names(by_chr)) {
    d <- copy(by_chr[[chr_name]])
    d[, idx := seq_len(.N)]
    runs <- d[, .(idx_start = min(idx), run_len = .N), by = run_id][order(idx_start)]
    out[[chr_name]] <- list(
      stat = d$stat,
      run_start = runs$idx_start,
      run_len = runs$run_len
    )
  }
  out
}

candidate_starts <- function(run_start, run_len, k) {
  ok <- which(run_len >= k)
  if (!length(ok)) return(integer())
  unlist(
    Map(
      function(s, l) seq.int(s, s + l - k),
      run_start[ok],
      run_len[ok]
    ),
    use.names = FALSE
  )
}

build_bg_candidate_medians <- function(bg_struct, regions) {
  out <- list()
  if (!length(bg_struct)) return(out)
  for (chr_name in names(bg_struct)) {
    if (!nrow(regions[chr == chr_name])) next
    needed_k <- sort(unique(regions[chr == chr_name, n_windows]))
    chr_stat <- bg_struct[[chr_name]]$stat
    chr_run_start <- bg_struct[[chr_name]]$run_start
    chr_run_len <- bg_struct[[chr_name]]$run_len
    out_chr <- list()
    for (k in needed_k) {
      starts <- candidate_starts(chr_run_start, chr_run_len, k)
      if (!length(starts)) {
        out_chr[[as.character(k)]] <- numeric()
      } else {
        meds <- vapply(
          starts,
          function(s) median(chr_stat[s:(s + k - 1L)], na.rm = TRUE),
          numeric(1)
        )
        out_chr[[as.character(k)]] <- meds
      }
    }
    out[[chr_name]] <- out_chr
  }
  out
}

run_one_combo <- function(dt, metric_name, pop_name, n_perm) {
  regions <- make_peak_regions(dt)
  bg_struct <- make_bg_struct(dt)
  if (!nrow(regions) || !length(bg_struct)) {
    return(list(
      summary = data.table(
        metric = metric_name,
        population = pop_name,
        n_peak_regions = 0L,
        n_perm = n_perm,
        n_perm_valid = 0L,
        obs_stat = NA_real_,
        null_median = NA_real_,
        delta_obs_minus_null = NA_real_,
        p_left = NA_real_,
        p_right = NA_real_,
        p_two_sided = NA_real_,
        min_matched_regions_per_perm = NA_integer_,
        median_matched_regions_per_perm = NA_real_
      ),
      regions = data.table(),
      perm = data.table()
    ))
  }

  bg_candidates <- build_bg_candidate_medians(bg_struct, regions)
  obs_stat <- median(regions$region_stat, na.rm = TRUE)

  perm_stats <- rep(NA_real_, n_perm)
  matched_n <- integer(n_perm)
  chr_vec <- regions$chr
  k_vec <- regions$n_windows
  n_regions <- nrow(regions)

  for (b in seq_len(n_perm)) {
    sampled <- rep(NA_real_, n_regions)
    for (i in seq_len(n_regions)) {
      chr_i <- chr_vec[i]
      k_i <- as.character(k_vec[i])
      cand <- bg_candidates[[chr_i]][[k_i]]
      if (length(cand)) {
        sampled[i] <- sample(cand, size = 1L)
      }
    }
    matched_n[b] <- sum(is.finite(sampled))
    if (matched_n[b] > 0L) {
      perm_stats[b] <- median(sampled, na.rm = TRUE)
    }
  }

  valid <- is.finite(perm_stats)
  n_valid <- sum(valid)
  null_median <- if (n_valid > 0L) median(perm_stats[valid]) else NA_real_

  if (n_valid > 0L && is.finite(obs_stat)) {
    p_left <- (1 + sum(perm_stats[valid] <= obs_stat)) / (n_valid + 1)
    p_right <- (1 + sum(perm_stats[valid] >= obs_stat)) / (n_valid + 1)
    p_two <- (1 + sum(abs(perm_stats[valid] - null_median) >= abs(obs_stat - null_median))) / (n_valid + 1)
  } else {
    p_left <- NA_real_
    p_right <- NA_real_
    p_two <- NA_real_
  }

  summary_dt <- data.table(
    metric = metric_name,
    population = pop_name,
    n_peak_regions = n_regions,
    n_perm = n_perm,
    n_perm_valid = n_valid,
    obs_stat = obs_stat,
    null_median = null_median,
    delta_obs_minus_null = obs_stat - null_median,
    p_left = p_left,
    p_right = p_right,
    p_two_sided = p_two,
    min_matched_regions_per_perm = if (length(matched_n)) min(matched_n) else NA_integer_,
    median_matched_regions_per_perm = if (length(matched_n)) median(matched_n) else NA_real_
  )

  region_dt <- copy(regions)[, `:=`(
    metric = metric_name,
    population = pop_name
  )]
  perm_dt <- data.table(
    metric = metric_name,
    population = pop_name,
    perm_id = seq_len(n_perm),
    perm_stat = perm_stats,
    matched_regions = matched_n,
    obs_stat = obs_stat
  )

  list(summary = summary_dt, regions = region_dt, perm = perm_dt)
}

message("Reading input files...")
tajd_all <- rbindlist(list(
  read_tajd(file.path(IN_DIR, "wild_labeled_tajD.txt"), "Wild"),
  read_tajd(file.path(IN_DIR, "land_labeled_tajD.txt"), "Landrace"),
  read_tajd(file.path(IN_DIR, "imp_labeled_tajD.txt"), "Improved")
), use.names = TRUE)

pi_all <- rbindlist(list(
  read_pi(file.path(IN_DIR, "wild_labeled_pi_withN.txt"), "Wild"),
  read_pi(file.path(IN_DIR, "land_labeled_pi_withN.txt"), "Landrace"),
  read_pi(file.path(IN_DIR, "imp_labeled_pi_withN.txt"), "Improved")
), use.names = TRUE)

all_dt <- rbindlist(list(tajd_all, pi_all), use.names = TRUE, fill = TRUE)
setorder(all_dt, metric, pop, chr, start, end)

combos <- unique(all_dt[, .(metric, pop)])
all_summary <- list()
all_regions <- list()
all_perm <- list()

for (i in seq_len(nrow(combos))) {
  m <- combos$metric[i]
  p <- combos$pop[i]
  message(sprintf("Running %s - %s (%d permutations)...", m, p, N_PERM))
  dt_sub <- all_dt[metric == m & pop == p]
  res <- run_one_combo(dt_sub, m, p, N_PERM)
  all_summary[[paste(m, p, sep = "_")]] <- res$summary
  all_regions[[paste(m, p, sep = "_")]] <- res$regions
  all_perm[[paste(m, p, sep = "_")]] <- res$perm
}

summary_dt <- rbindlist(all_summary, fill = TRUE)
regions_dt <- rbindlist(all_regions, fill = TRUE)
perm_dt <- rbindlist(all_perm, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "matched_permutation_summary.csv"))
fwrite(regions_dt, file.path(OUT_DIR, "matched_peak_regions_observed.csv"))
fwrite(perm_dt, file.path(OUT_DIR, "matched_permutation_null_stats.csv"))

plot_dt <- merge(
  perm_dt[is.finite(perm_stat)],
  summary_dt[, .(metric, population, p_two_sided)],
  by = c("metric", "population"),
  all.x = TRUE
)

plot_dt[, panel := sprintf("%s - %s", metric, population)]
plot_dt[, panel := factor(
  panel,
  levels = c(
    "TajimaD - Wild", "TajimaD - Landrace", "TajimaD - Improved",
    "Pi - Wild", "Pi - Landrace", "Pi - Improved"
  )
)]

obs_lines <- unique(plot_dt[, .(metric, population, panel, obs_stat)])
ann <- unique(plot_dt[, .(metric, population, panel, p_two_sided)])
ann[, label := sprintf("p(two-sided)=%.4g", p_two_sided)]

p <- ggplot(plot_dt, aes(x = perm_stat)) +
  geom_histogram(bins = 60, fill = "#BFBFBF", color = "white") +
  geom_vline(
    data = obs_lines,
    aes(xintercept = obs_stat),
    color = "#B22222",
    linewidth = 0.8
  ) +
  geom_text(
    data = ann,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.02, vjust = 1.5, size = 3
  ) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  labs(
    x = "Permutation null statistic (median of matched-region medians)",
    y = "Count",
    title = sprintf("Matched permutation test (N=%d)", N_PERM),
    subtitle = "Red line = observed statistic from real PEAK regions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUT_DIR, "matched_permutation_null_histograms.pdf"),
  p, width = 13, height = 8.5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "matched_permutation_null_histograms.png"),
  p, width = 13, height = 8.5, dpi = 300
)

message("Done.")
message("Output directory: ", OUT_DIR)
