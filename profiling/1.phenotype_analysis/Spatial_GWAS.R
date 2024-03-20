library(SpATS)
library(tidyverse)
library(rMVP)

trait="PH"

phenotype <- read.csv("data/Revised_exp_data/DW_Revised.csv")
#phenotype[, 6:8] <- sapply(phenotype[, 6:8], as.numeric)
mymean <- function(x) mean(x, na.rm=T);
phenotype$Dry_shoot_weight <- apply(phenotype[,6:8], 1, mymean)


#phenotype$plant_height_rep_3 <- as.numeric(phenotype$plant_height_rep_3)
#phenotype <- pivot_longer(phenotype, c(plant_height_rep_1, plant_height_rep_2,plant_height_rep_3), 
                          #names_to = "rep", values_to = trait)

phenotype$Treatment <- as.factor(phenotype$Treatment)
phenotype$rep <- as.factor(phenotype$rep)
phenotype$C <- as.factor(phenotype$Column)
phenotype$R <- as.factor(phenotype$Row)
phenotype$Pedigree <- as.factor(phenotype$Pedigree)
phenotype$Pot_id <- as.factor(phenotype$Pot_id)

#phenotype_avg <- phenotype %>% group_by(across(everything())) %>% summarise(PlantHeight = mean(PH, na.rm = TRUE)) %>%  ungroup()

m0 <- SpATS(response = "Dry_shoot_weight", spatial = ~ SAP(Column, Row, nseg = c(9,20)),
            genotype = "Pedigree", fixed = ~ Treatment , random = ~ R + C +Pot_id,
            #genotype.as.random = TRUE,
            data = phenotype, control = list(tolerance = 1e-03))


m_avg <- SpATS(response = "Dry_shoot_weight", spatial = ~ SAP(Column, Row, nseg = c(9,20)),
            genotype = "Pedigree", fixed = ~ Treatment, random = ~ R + C ,
            #genotype.as.random = TRUE,
            data = phenotype, control = list(tolerance = 1e-03))
plot(m_avg)


phenotype_hn <- phenotype[phenotype$Treatment=="HN",]
phenotype_ln <- phenotype[phenotype$Treatment=="LN",]
m_hn <- SpATS(response = "mval", spatial = ~ SAP(Column, Row, nseg = c(9,20)),
               genotype = "Pedigree", random = ~ R + C ,
               #genotype.as.random = TRUE,
               data = phenotype_hn, control = list(tolerance = 1e-03))
m_ln <- SpATS(response = "mval", spatial = ~ SAP(Column, Row, nseg = c(9,20)),
               genotype = "Pedigree", random = ~ R + C ,
               #genotype.as.random = TRUE,
               data = phenotype_ln, control = list(tolerance = 1e-03))


# Open a PNG device
png(filename = "graphs/saptial graphs/m_hn_plot.png", width = 800, height = 600, units = "px", res = 500)
plot(m_hn) # Assuming m_hn is the object you want to plot
dev.off()
# Saving the second plot as PNG
png("path/to/save/m_ln_plot.png", width=800, height=600)
plot(m_ln) # Assuming m_ln is the data/object you want to plot
dev.off() # Close the device


blue_hn <- as.data.frame(m_hn$coeff[1:343]+m_hn$coeff[344])
blue_ln <- as.data.frame(m_ln$coeff[1:345]+m_ln$coeff[346])

write.table(blue_hn,paste0("data/BLUES/DW_hn_blue.csv"), sep=",", quote=F)
write.table(blue_ln,paste0("data/BLUES/DW_ln_blue.csv"), sep=",", quote=F)
