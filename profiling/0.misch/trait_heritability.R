library(lme4)
library(tidyr)
library(dplyr)
library(ggplot2)

df <- read.csv("1.2024_SAP_Pheno_raw.csv")
ph_long <- df %>%
  select(geno, block, Plant.Height.1, Plant.Height.2, Plant.Height.3) %>%
  pivot_longer(cols = starts_with("Plant.Height"),
               names_to = "replicate",
               values_to = "plant_height")
str(df)
head(df)
head(ph_long)
# Make sure geno and replicate are factors
ph_long$geno <- as.factor(ph_long$geno)
ph_long$replicate <- as.factor(ph_long$replicate)
# Fit model with genotype as random effect
model <- lmer(plant_height ~ (1|geno) + (1|replicate), data = ph_long)
# Extract variance components
varcomp <- as.data.frame(VarCorr(model))
vg <- varcomp$vcov[varcomp$grp == "geno"]         # genetic variance
ve <- sigma(model)^2                              # residual variance
n_rep <- 3                                         # 3 replicates
# Broad-sense heritability
H2 <- vg / (vg + ve / n_rep)
print(H2)




chl_long <- df %>%
  select(geno, block, Chl.1, Chl.2, Chl.3) %>%
  pivot_longer(cols = starts_with("Chl."),
               names_to = "replicate",
               values_to = "chlorophyll")

chl_long$geno <- as.factor(chl_long$geno)
chl_long$replicate <- as.factor(chl_long$replicate)
chl_model <- lmer(chlorophyll ~ (1|geno) + (1|replicate), data = chl_long)
# Extract variance components
varcomp_chl <- as.data.frame(VarCorr(chl_model))
vg_chl <- varcomp_chl$vcov[varcomp_chl$grp == "geno"]
ve_chl <- sigma(chl_model)^2
n_rep <- 3
H2_chl <- vg_chl / (vg_chl + ve_chl / n_rep)
print(H2_chl)


ph_long_labeled <- ph_long %>%
  rename(value = plant_height) %>%
  mutate(trait = "Plant Height")

chl_long_labeled <- chl_long %>%
  rename(value = chlorophyll) %>%
  mutate(trait = "Chlorophyll")


combined_traits <- bind_rows(ph_long_labeled, chl_long_labeled)


ggplot(combined_traits, aes(x = trait, y = value, fill = trait)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  labs(title = "Heritability comparison of plant height and chlorophyll",
       x = "Trait", y = "Heritability") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("/Users/subhashmahamkali/Documents/gwas_sap/data/1.Phenotype_data/3.2024/trait_heritability_comparison.pdf",
       width = 6, height = 5, units = "in",
       #device = cairo_pdf,  # better text rendering
       bg = "white") 
