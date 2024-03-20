#https://speciationgenomics.github.io/pca/
#pca was conducted using VCF file susing plink on HCC.

library(tidyverse)
# read in data
pca <- read_table2("data/sorghum.eigenvec", col_names = FALSE)
eigenval <- scan("data/sorghum.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "id"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
library(ggplot2)
ggplot(data=pca, aes(PC1,PC2)) + geom_point()

