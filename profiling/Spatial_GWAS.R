library(SpATS)
library(tidyverse)
library(rMVP)

trait="PH"

phenotype <- read.csv(paste0(trait,".csv"))
phenotype$plant_height_rep_3 <- as.numeric(phenotype$plant_height_rep_3)
phenotype <- pivot_longer(phenotype, c(plant_height_rep_1, plant_height_rep_2,plant_height_rep_3), 
                          names_to = "rep", values_to = trait)
phenotype$Treatment <- as.factor(phenotype$Treatment)
phenotype$rep <- as.factor(phenotype$rep)
phenotype$C <- as.factor(phenotype$Column)
phenotype$R <- as.factor(phenotype$Row)
phenotype$Pedigree <- as.factor(phenotype$Pedigree)
phenotype$Pot_id <- as.factor(phenotype$Pot_id)

phenotype_avg <- phenotype %>% group_by(across(everything())) %>% ungroup(PlantHeight, rep) %>%
  summarise(PlantHeight=mean(PlantHeight, na.rm=T))

m0 <- SpATS(response = trait, spatial = ~ SAP(Column, Row, nseg = c(9,20)),
            genotype = "Pedigree", fixed = ~ Treatment , random = ~ R + C +Pot_id,
            #genotype.as.random = TRUE,
            data = phenotype, control = list(tolerance = 1e-03))

m_avg <- SpATS(response = trait, spatial = ~ SAP(Column, Row, nseg = c(9,20)),
            genotype = "Pedigree", fixed = ~ Treatment, random = ~ R + C ,
            #genotype.as.random = TRUE,
            data = phenotype_avg, control = list(tolerance = 1e-03))
plot(m_avg)

phenotype_hn <- phenotype_avg[phenotype_avg$Treatment=="HN",]
phenotype_ln <- phenotype_avg[phenotype_avg$Treatment=="LN",]
m_hn <- SpATS(response = trait, spatial = ~ SAP(Column, Row, nseg = c(9,20)),
               genotype = "Pedigree", random = ~ R + C ,
               #genotype.as.random = TRUE,
               data = phenotype_hn, control = list(tolerance = 1e-03))
m_ln <- SpATS(response = trait, spatial = ~ SAP(Column, Row, nseg = c(9,20)),
               genotype = "Pedigree", random = ~ R + C ,
               #genotype.as.random = TRUE,
               data = phenotype_ln, control = list(tolerance = 1e-03))
plot(m_hn)#343 samples
plot(m_ln)#345 samples

blue_hn <- as.data.frame(m_hn$coeff[1:343]+m_hn$coeff[344])
blue_ln <- as.data.frame(m_ln$coeff[1:345]+m_ln$coeff[346])

write.table(blue_hn,paste0(trait, "_hn_blue.csv"), sep=",", quote=F)
write.table(blue_ln,paste0(trait, "_ln_blue.csv"), sep=",", quote=F)
