library(data.table)
library(dplyr)

# Load data
individual_data <- fread("/work/jyanglab/subhash/sorgsd/filterfed_vcf/genotypes_fixed.txt", sep = "\t")
individual_data$V1 <- ifelse(individual_data$V1 %in% 1:9, paste0("Chr0", individual_data$V1), paste0("Chr", individual_data$V1))
setnames(individual_data, c("chr", "pos", "ref", "alt", paste0("indiv", 1:(ncol(individual_data) - 4))))

overlapping_data <- fread("/work/jyanglab/subhash/sorgsd/filterfed_vcf/0.vcf_chunk/0.unziped/tsv_out/0.zero_shot/1.bed_files/zeroshot_in_gene.bed")
setnames(overlapping_data, c("chr", "start", "end", "pos", "zeroshotscore", "chr_1", "start_1", "end_1", "geneID", "annot", "string", "int"))

# Merge
merged_data <- merge(overlapping_data, individual_data, by = c("chr", "pos"))
d <- merged_data[, c("chr", "pos", "zeroshotscore", "start_1", "end_1", "geneID", paste0("indiv", 1:289))]

# Filter for deleterious SNPs
subset_df <- d[zeroshotscore < 0, ]

# --- VECTORIZE GENETIC LOAD ---

# Function to convert genotype to multiplier
geno_to_multiplier <- function(x) {
  ifelse(x == "./.", NA,
    ifelse(x == "1/1", 2,
      ifelse(x %in% c("0/1", "1/0"), 1,
        ifelse(x == "0/0", 0, NA)
      )
    )
  )
}

# Apply to all genotype columns at once
genotype_cols <- paste0("indiv", 1:289)
genotype_matrix <- subset_df[, ..genotype_cols]

# Convert all genotypes to multipliers
multiplier_matrix <- as.data.table(lapply(genotype_matrix, geno_to_multiplier))

# Multiply by zeroshotscore
for (i in seq_along(multiplier_matrix)) {
  multiplier_matrix[[i]] <- multiplier_matrix[[i]] * subset_df$zeroshotscore
}

# Combine gene info + load info
subset_df_final <- cbind(subset_df[, .(chr, pos, start_1, end_1, geneID)], multiplier_matrix)

# --- SUMMARIZE ---
# Create genetic load columns names
genetic_load_cols <- names(multiplier_matrix)

# Summarize per gene
gene_load_summary <- subset_df_final[, c(
  .(chr = first(chr), start = min(start_1), end = max(end_1)),
  lapply(.SD, function(x) sum(abs(x), na.rm = TRUE))
), by = geneID, .SDcols = genetic_load_cols]

# Save output
fwrite(gene_load_summary, "/work/jyanglab/subhash/sorgsd/filterfed_vcf/gene_load_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
