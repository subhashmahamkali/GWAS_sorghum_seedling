library(data.table)
library(readxl)
mor = read_excel("/Users/subhashmahamkali/Downloads/sd01 (1).xlsx")
mor = mor[, 1]
mor <- mor[-1, ]
colnames(mor)[1] <- "sample"
gsd = read_excel("/Users/subhashmahamkali/Downloads/13068_2021_2016_MOESM1_ESM (2).xlsx")
gsd = gsd[, 1:2]
gsd = gsd[-1, ]
gsd[[2]] <- gsub(" ", "", gsd[[2]])
colnames(gsd)[1] = "sample"
colnames(gsd)[2] = "plant"
sap = fread("/Users/subhashmahamkali/Downloads/tpj15853-sup-0001-files1 (2).csv")
sap = sap[, 1]
colnames(sap)[1] = "sample"
mor_sap = intersect(mor$sample, sap$sample)
mor_sap = data.frame(sample = mor_sap)
mor_gsd = intersect(mor$sample, gsd$plant)
mor_gsd = data.frame(sample = mor_gsd)



d = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2_gene.bed")
b <- read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/7.sorghum.annotations.xlsx", skip = 2)
b <- b[, c(1, 9)]
b <- b[-c(1, 2), ]
b <- b[!duplicated(b$Gene), ]
colnames(b) <- c("Gene", "Description")
colnames(d)[9] <- "Gene"
merged_data <- merge(d, b, by = "Gene", all.x = TRUE)
write.table(merged_data, "/Users/subhashmahamkali/Downloads/1.miscellaneous/b2_genes_annot.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)




library(data.table)
# Load original score intervals
orig <- fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2.bed", header = FALSE)
colnames(orig) <- c("chr", "start", "end", "score")
# Load gene-intersected file
gene <- fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2_gene.bed", header = FALSE)
colnames(gene) <- c("chr", "start", "end", "score", "gene_chr", "gene_start", "gene_end", "strand", "gene_id")

# For each gene-intersected row, find the matching original score region (exact max score)
gene_with_orig_coords <- gene[, {
  # Filter orig to match the current row's chr and overlapping coordinates
  matching_orig <- orig[chr == .SD$chr & end >= .SD$start & start <= .SD$end]
  
  # Find the one with score closest to the gene's intersect score
  hit <- matching_orig[which.min(abs(score - .SD$score))]
  
  if (nrow(hit) == 1 && abs(hit$score - .SD$score) < 1e-6) {
    .SD[, .(chr, start, end, score, gene_chr, gene_start, gene_end, strand, gene_id,
            orig_start = hit$start, orig_end = hit$end)]
  } else {
    NULL
  }
}, by = 1:nrow(gene)]  # use row index for safe row-wise ops

# Save the result
fwrite(gene_with_orig_coords,
       file = "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2_gene_with_orig_coords.tsv",
       sep = "\t", quote = FALSE)


library(data.table)
# Step 1: Load merged gene file (with max_score per merged region)
gene <- fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/pos_genes.bed", header = FALSE)
colnames(gene) <- c("chr", "merged_start", "merged_end", "max_score", "comparison",
                    "gene_chr", "gene_start", "gene_end", "strand", "gene_id")

# Step 2: Load original FST files into a named list
fst_files <- list(
  sorg_wild_vs_improved = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg_wild_vs_improved.windowed.weir.fst"),
  landrace_vs_improved  = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/landrace_vs_improved.windowed.weir.fst"),
  sorg_wild_vs_landrace = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/sorg_wild_vs_landrace.windowed.weir.fst")
)
# Ensure proper column names
for (name in names(fst_files)) {
  setnames(fst_files[[name]], c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST"))
}

# Step 3: Apply fuzzy match per comparison
result_list <- list()
for (comp in unique(gene$comparison)) {
  gene_sub <- gene[comparison == comp]
  fst <- fst_files[[comp]]
  result_sub <- gene_sub[, {
    hits <- fst[CHROM == chr & BIN_END >= merged_start & BIN_START <= merged_end]
    match <- hits[which.min(abs(WEIGHTED_FST - max_score))]
    if (nrow(match) == 1) {
      .(chr, merged_start, merged_end, max_score, comparison,
        gene_chr, gene_start, gene_end, strand, gene_id,
        fst_exact_start = match$BIN_START,
        fst_exact_end   = match$BIN_END)
    } else {
      NULL
    }
  }, by = 1:nrow(gene_sub)]
  
  result_list[[comp]] <- result_sub
}

# Step 4: Combine all results
final_result <- rbindlist(result_list, use.names = TRUE)

# Step 5: Save to file
fwrite(final_result,
       "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/merged_genes_with_fst_exact_coords_all_comparisons.tsv",
       sep = "\t", quote = FALSE)



b <- read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/7.sorghum.annotations.xlsx", skip = 2)
b <- b[, c(1, 9)]
b <- b[-c(1, 2), ]
b <- b[!duplicated(b$Gene), ]

colnames(b) <- c("Gene", "Description")
#d = gene_with_orig_coords
#colnames(d)[10] <- "Gene"
merged_data <- merge(d, b, by = "Gene", all.x = TRUE)
write.table(merged_data, "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/B2/b2_genes_annot.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#d = final_result
#colnames(d)[11] <- "Gene"
write.table(merged_data, "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/positive_selection/pos_genes_annot.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)



#after this I merged these 2 files on excell with b2 as column and all these as another column and saved it.



d = fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum_project/GWAS/GWAS_NR.bed")
max(d$V4)
min(d$V4)

