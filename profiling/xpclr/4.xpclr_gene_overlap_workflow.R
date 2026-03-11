library(data.table)

# XP-CLR gene-overlap workflow
# Starting inputs are the merged genome-wide XP-CLR result files produced on HCC.

base_dir <- "/work/jyanglab/subhash/sorgsd/xpclr"
out_dir <- file.path(base_dir, "processed_candidates")
gene_bed <- "/mnt/nrdstor/jyanglab/subhash/datasets/3.results/6.BT623_gene.V5_up_dw_2k.bed"
merge_dist <- 40000L

dir.create(file.path(out_dir, "01_all_windows_bed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "02_top1pct"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "03_sorted"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "04_merged"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "05_gene_overlap"), recursive = TRUE, showWarnings = FALSE)

comparisons <- list(
  domestication_landrace_vs_wild = file.path(
    base_dir,
    "domestication_landrace_vs_wild",
    "domestication_landrace_vs_wild.all_chr.tsv"
  ),
  breeding_improved_vs_landrace = file.path(
    base_dir,
    "breeding_improved_vs_landrace",
    "breeding_improved_vs_landrace.all_chr.tsv"
  )
)

make_xpclr_bed <- function(infile, label) {
  d <- fread(infile)

  required_cols <- c(
    "chrom", "start", "stop", "pos_start", "pos_stop",
    "nSNPs", "nSNPs_avail", "xpclr", "xpclr_norm"
  )
  stopifnot(all(required_cols %in% names(d)))

  d[, chrom := as.integer(chrom)]
  d[, start := as.integer(start)]
  d[, stop := as.integer(stop)]
  d[, pos_start := as.integer(pos_start)]
  d[, pos_stop := as.integer(pos_stop)]
  d[, nSNPs := as.integer(nSNPs)]
  d[, nSNPs_avail := as.integer(nSNPs_avail)]
  d[, xpclr := as.numeric(xpclr)]
  d[, xpclr_norm := as.numeric(xpclr_norm)]

  bed <- d[, .(
    chr = paste0("Chr", sprintf("%02d", chrom)),
    start = start,
    end = stop,
    pos_start = pos_start,
    pos_stop = pos_stop,
    nSNPs = nSNPs,
    nSNPs_avail = nSNPs_avail,
    xpclr = xpclr,
    xpclr_norm = xpclr_norm
  )]

  all_out <- file.path(out_dir, "01_all_windows_bed", paste0(label, ".all_windows.bed"))
  fwrite(bed, all_out, sep = "\t", col.names = FALSE)

  thr <- quantile(bed$xpclr, 0.99, na.rm = TRUE)
  top <- bed[xpclr >= thr]

  top_out <- file.path(out_dir, "02_top1pct", paste0(label, ".top1pct.bed"))
  thr_out <- file.path(out_dir, "02_top1pct", paste0(label, ".top1pct.threshold.txt"))

  fwrite(top, top_out, sep = "\t", col.names = FALSE)
  fwrite(
    data.table(comparison = label, top1pct_threshold = thr, n_windows = nrow(top)),
    thr_out,
    sep = "\t"
  )

  invisible(list(top_out = top_out, threshold = thr))
}

sort_top_bed <- function(label) {
  top_file <- file.path(out_dir, "02_top1pct", paste0(label, ".top1pct.bed"))
  sorted_file <- file.path(out_dir, "03_sorted", paste0(label, ".top1pct.sorted.bed"))

  top <- fread(top_file, header = FALSE)
  setorder(top, V1, V2, V3)
  fwrite(top, sorted_file, sep = "\t", col.names = FALSE)

  invisible(sorted_file)
}

merge_top_regions <- function(label) {
  sorted_file <- file.path(out_dir, "03_sorted", paste0(label, ".top1pct.sorted.bed"))
  merged_file <- file.path(out_dir, "04_merged", paste0(label, ".top1pct.merged_max_count.bed"))

  top <- fread(sorted_file, header = FALSE)
  setnames(top, c("chr", "start", "end", "pos_start", "pos_stop", "nSNPs", "nSNPs_avail", "xpclr", "xpclr_norm"))

  result <- vector("list", 0)
  i <- 1L
  n <- nrow(top)

  while (i <= n) {
    chr <- top$chr[i]
    region_start <- top$start[i]
    region_end <- top$end[i]
    max_xpclr <- top$xpclr[i]
    count_windows <- 1L
    j <- i + 1L

    while (j <= n && identical(top$chr[j], chr) && top$start[j] - region_end <= merge_dist) {
      region_end <- max(region_end, top$end[j])
      max_xpclr <- max(max_xpclr, top$xpclr[j], na.rm = TRUE)
      count_windows <- count_windows + 1L
      j <- j + 1L
    }

    result[[length(result) + 1L]] <- data.table(
      chr = chr,
      merged_start = region_start,
      merged_end = region_end,
      max_xpclr = max_xpclr,
      n_merged_windows = count_windows
    )

    i <- j
  }

  merged <- rbindlist(result)
  fwrite(merged, merged_file, sep = "\t", col.names = FALSE)

  invisible(merged_file)
}

fix_chr_names <- function(label) {
  merged_file <- file.path(out_dir, "04_merged", paste0(label, ".top1pct.merged_max_count.bed"))
  chrfix_file <- file.path(out_dir, "04_merged", paste0(label, ".top1pct.merged_max_count.chrfix.bed"))

  x <- fread(merged_file, header = FALSE)
  x[, V1 := sub("^Chr0", "", V1)]
  x[, V1 := sub("^Chr", "", V1)]
  fwrite(x, chrfix_file, sep = "\t", col.names = FALSE)

  invisible(chrfix_file)
}

overlap_with_genes <- function(label) {
  chrfix_file <- file.path(out_dir, "04_merged", paste0(label, ".top1pct.merged_max_count.chrfix.bed"))
  overlap_bed <- file.path(out_dir, "05_gene_overlap", paste0(label, ".top1pct.merged_max_count_with_genes.bed"))
  overlap_txt <- file.path(out_dir, "05_gene_overlap", paste0(label, ".top1pct.merged_max_count_with_genes.txt"))

  a <- fread(chrfix_file, header = FALSE)
  b <- fread(gene_bed, header = FALSE)

  setnames(a, c("region_chr", "region_start", "region_end", "max_xpclr", "n_merged_windows"))
  setnames(b, c("gene_chr", "gene_start", "gene_end", "strand", "gene_id"))

  a[, join_chr := as.character(region_chr)]
  b[, join_chr := as.character(gene_chr)]

  overlap <- a[b, on = .(join_chr)]
  overlap <- overlap[
    !(region_end < gene_start | region_start > gene_end),
    .(
      region_chr,
      region_start,
      region_end,
      max_xpclr,
      n_merged_windows,
      gene_chr,
      gene_start,
      gene_end,
      strand,
      gene_id
    )
  ]

  setorder(overlap, region_chr, region_start, region_end, gene_start, gene_end)
  overlap <- unique(overlap)

  fwrite(overlap, overlap_bed, sep = "\t", col.names = FALSE)

  overlap[, comparison := label]
  setcolorder(
    overlap,
    c(
      "comparison", "region_chr", "region_start", "region_end",
      "max_xpclr", "n_merged_windows", "gene_chr", "gene_start",
      "gene_end", "strand", "gene_id"
    )
  )
  fwrite(overlap, overlap_txt, sep = "\t")

  invisible(list(bed = overlap_bed, txt = overlap_txt))
}

for (label in names(comparisons)) {
  message("Processing: ", label)
  make_xpclr_bed(comparisons[[label]], label)
  sort_top_bed(label)
  merge_top_regions(label)
  fix_chr_names(label)
  overlap_with_genes(label)
}
