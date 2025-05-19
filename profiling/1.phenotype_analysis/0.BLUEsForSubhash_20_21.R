#Author: jensina


library(tidyverse)
library(SpATS)
library(readxl)
library(lme4)

# Functions
printHistogram <- function(data, phenotype, title = NULL, bin_no=50)
{
  df <- data
  p <- ggplot(df, aes(.data[[phenotype]])) +
    geom_histogram(bins = bin_no) +
    labs(title = title)
  print(p)
}

# Use new, fixed version of data from James
# Rename columns because I'm picky (and it makes combining across year datasets easier later)
# Fix genotypes with multiple names in the dataset 
sap2020 <- read_excel('/Users/subhashmahamkali/Downloads/subhash-blues/SAP2020_merged_v3.1.xls') %>%
  rename(plot = PlotID, 
         genotype = SorghumAccession,
         poorStand = `PoorStand?`, 
         row = Row, 
         column = Column, 
         block = Block, 
         daysToFlower = DaysToBloom, 
         medianLeafAngle = MedianLeafAngle,
         leafAngleStandardDeviation = LeafAngleSDV,
         paniclesPerPlot = PaniclesPerPlot, 
         panicleGrainWeight = PanicleGrainWeight,
         estimatedPlotYield = EstimatedPlotYield, 
         flagLeafLength = FlagLeafLength, 
         flagLeafWidth = FlagLeafWidth,
         extantLeafNumber = ExtantLeafNumber, 
         plantHeight = PlantHeight, 
         thirdLeafLength = ThirdLeafLength,
         thirdLeafWidth = ThirdLeafWidth,
         tillersPerPlant = TillersPerPlant,
         stemDiameterLower = StemDiameterLower,
         stemDiameterUpper = StemDiameterUpper,
         rachisLength = RachisLength, 
         rachisDiameterLower = RachisDiameterLower,
         rachisDiameterUpper = RachisDiameterUpper,
         primaryBranchNumber = PrimaryBranchNo, 
         branchInternodeLength = BranchInternodeLength,
         percentMoisture = MoisturePCT,
         percentProtein = ProteinPCT,
         percentOil = OilPCT,
         percentAsh = AshPCT, 
         percentStarch = StarchPCT, 
         seedColor = KernelColor) %>%
  rowwise() %>%
  mutate(genotype = case_when(genotype=="Btx623"|genotype=="BTx623" ~ 'PI 564163', 
                              genotype=='TX430' ~ 'PI 655996', 
                              .default = genotype), 
         treatment = case_when(Treatment=='LowNitrogen' ~ 'LN', 
                               Treatment=='SufficientNitrogen' ~ 'HN'),
         plotLocation = str_c(row, column, block, sep = '-')) %>%
  select(!c(SorghumName, SNPDataID, Treatment))
sap2020kyle <- read.csv('/Users/subhashmahamkali/Downloads/subhash-blues/2020_Sorghum_KL_Data.csv') %>%
  rowwise() %>%
  mutate(plotLocation = str_c(Row, Col, Rep, sep='-'))

sap2020 <- filter(sap2020, block==4|plotLocation %in% sap2020kyle$plotLocation) %>%
  filter(!is.na(column))
#Split dataset into high and low N
sap2020LN <- sap2020 %>% 
  filter(treatment=='LN')
sap2020HN <- sap2020 %>% 
  filter(treatment=='HN')

# Repeat process for 2021 (aside from Kyle's corrections)
sap2021 <- read_excel('/Users/subhashmahamkali/Downloads/subhash-blues/SAPMerged2021_v2.3.xlsx')[1:27] %>% 
  rename(genotype = PINumber, 
         plot = Plot,
         treatment = Treatment, 
         block = Rep, 
         row = Row, 
         column = Column, 
         plantHeight = PlantHeight, 
         daysToFlower = DaysToFlower, 
         thirdLeafWidth = LeafWidth, 
         stemDiameterLower = StemDiameterLower, 
         stemDiameterUpper = StemDiameterUpper,
         paniclesPerPlant = PaniclesPerPlant, 
         extantLeafNumber = LeafNumber, 
         thirdLeafLength = LeafLength, 
         seedMassPerPlant = SeedMassPerPlant, 
         percentProtein = SeedProtein, 
         percentAsh = SeedAsh,
         percentOil = SeedOil,
         percentStarch = SeedStarch, 
         percentMoisture = SeedMoisture, 
         seedColor = SeedColor, 
         tillersPerPlant = TillersPerPlant) %>%
  mutate(plantHeight = plantHeight*100) %>%
  rowwise() %>%
  mutate(genotype = case_when(genotype=="Tx430"|genotype=='RTx430' ~ 'PI 655996',
                              .default = genotype)) %>%
  filter(genotype %in% sap2020$genotype) %>%
  select(!c(GenoIDFromMiaoEtAl, SAPID, SorghumConversionID, Name, Comments))

sap2021LN <- sap2021 %>%
  filter(treatment=='LN')

sap2021HN <- sap2021 %>%
  filter(treatment=='HN')

envDatasets <- list(sap2020LN, sap2020HN, sap2021LN, sap2021HN)
envYears <- c('2020', '2020', '2021', '2021')
envNitrogen <- c('LN', 'HN', 'LN', 'HN')

phenotypes2020 <- colnames(sap2020)[c(6:8, 10:31)]
phenotypes2021 <- colnames(sap2021)[7:21]
phenotypeList <- list(phenotypes2020, phenotypes2020, phenotypes2021, phenotypes2021)

# Look at distributions to determine threshold of outliers to remove
for (i in 1:length(envDatasets))
{
  df <- envDatasets[i][[1]]
  phenotypes <- phenotypeList[i][[1]]
  for (j in phenotypes)
  {
    if(is.null(df[[j]]))
    {
      next
    }
    printHistogram(data = df, phenotype = j, title = paste0(envYears[i], ':', envNitrogen[i], ':', j))
  }
}
# Remove outliers 
sap2021HN <- sap2021HN %>%
  rowwise() %>%
  mutate(percentMoisture = case_when(percentMoisture > 14.5 ~ NA, .default = percentMoisture), 
         percentStarch = case_when(percentStarch < 55 ~ NA, .default = percentStarch), 
         percentAsh = case_when(percentAsh > 2.7 ~ NA, .default = percentAsh),
         percentOil = case_when(percentOil > 5.6 ~ NA, .default = percentOil),
         percentProtein = case_when(percentProtein < 5|percentProtein > 27 ~ NA, .default = percentProtein),
         seedMassPerPlant = case_when(seedMassPerPlant > 135 ~ NA, .default = seedMassPerPlant), 
         thirdLeafWidth = case_when(thirdLeafWidth > 20 ~ NA, .default = thirdLeafWidth),
         thirdLeafLength = case_when(thirdLeafLength < 28 | thirdLeafLength > 90 ~ NA, .default = thirdLeafLength),
         extantLeafNumber = case_when(extantLeafNumber < 5 | extantLeafNumber > 17 ~ NA, .default = extantLeafNumber),
         paniclesPerPlant = case_when(paniclesPerPlant > 4 ~ NA, .default = paniclesPerPlant), 
         stemDiameterUpper = case_when(stemDiameterUpper < 5 ~ NA, .default = stemDiameterUpper), 
         stemDiameterLower = case_when(stemDiameterLower < 10 | stemDiameterLower > 30 ~ NA, .default = stemDiameterLower), 
         tillersPerPlant = case_when(tillersPerPlant > 3 ~ NA, .default = tillersPerPlant), 
         plantHeight = case_when(plantHeight > 375 ~ NA, .default = plantHeight))

sap2021LN <- sap2021LN %>%
  rowwise() %>%
  mutate(stemDiameterUpper = case_when(stemDiameterUpper > 16 ~ NA, .default = stemDiameterUpper), 
         stemDiameterLower = case_when(stemDiameterLower > 31 ~ NA, .default = stemDiameterLower), 
         tillersPerPlant = case_when(tillersPerPlant > 2.5 ~ NA, .default = tillersPerPlant), 
         daysToFlower = case_when(daysToFlower < 57 | daysToFlower > 95 ~ NA, .default = daysToFlower), 
         plantHeight = case_when(plantHeight > 200 ~ NA, .default = plantHeight))

sap2020HN <- sap2020HN %>%
  rowwise() %>%
  mutate(percentStarch = case_when(percentStarch < 52 ~ NA, .default = percentStarch), 
         percentAsh = case_when(percentAsh > 2.15 ~ NA, .default = percentAsh),
         percentOil = case_when(percentOil > 6.25 ~ NA, .default = percentOil),
         percentProtein = case_when(percentProtein < 8 | percentProtein > 18 ~ NA, .default = percentProtein), 
         percentMoisture = case_when(percentMoisture < 6 | percentMoisture > 14 ~ NA, .default = percentMoisture),
         branchInternodeLength = case_when(branchInternodeLength > 10 ~ NA, .default = branchInternodeLength), 
         primaryBranchNumber = case_when(primaryBranchNumber > 130 ~ NA, .default = primaryBranchNumber),
         rachisDiameterUpper = case_when(rachisDiameterUpper > 6.25 ~ NA, .default = rachisDiameterUpper),
         rachisDiameterLower = case_when(rachisDiameterLower < 3 | rachisDiameterLower > 12 ~ NA, .default = rachisDiameterLower),
         rachisLength = case_when(rachisLength > 40 ~ NA, .default = rachisLength), 
         stemDiameterUpper = case_when(stemDiameterUpper < 5 | stemDiameterUpper > 15.5 ~ NA, .default = stemDiameterUpper),
         tillersPerPlant = case_when(tillersPerPlant > 10 ~ NA, .default = tillersPerPlant),
         thirdLeafWidth = case_when(thirdLeafWidth < 4 | thirdLeafWidth > 10 ~ NA, .default = thirdLeafWidth),  
         thirdLeafLength = case_when(thirdLeafLength < 40 ~ NA, .default = thirdLeafLength),
         plantHeight = case_when(plantHeight > 200 ~ NA, .default = plantHeight),
         flagLeafLength = case_when(flagLeafLength > 65 ~ NA, .default = flagLeafLength),
         estimatedPlotYield = case_when(estimatedPlotYield > 1500 ~ NA, .default = estimatedPlotYield),
         panicleGrainWeight = case_when(panicleGrainWeight > 125 ~ NA, .default = panicleGrainWeight),
         paniclesPerPlot = case_when(paniclesPerPlot > 50 ~ NA, .default = paniclesPerPlot),
         leafAngleStandardDeviation = case_when(leafAngleStandardDeviation > 15 ~ NA, .default = leafAngleStandardDeviation), 
         medianLeafAngle = case_when(medianLeafAngle < 20 ~ NA, .default = medianLeafAngle), 
         daysToFlower = case_when(daysToFlower > 80 ~ NA, .default = daysToFlower))

sap2020LN <- sap2020LN %>%
  rowwise() %>%
  mutate(branchInternodeLength = case_when(branchInternodeLength > 10 ~ NA, .default = branchInternodeLength), 
         primaryBranchNumber = case_when(primaryBranchNumber > 112 ~ NA, .default = primaryBranchNumber),
         rachisDiameterUpper = case_when(rachisDiameterUpper > 5 ~ NA, .default = rachisDiameterUpper),
         rachisDiameterLower = case_when(rachisDiameterLower > 12 ~ NA, .default = rachisDiameterLower),
         rachisLength = case_when(rachisLength < 3 | rachisLength > 40 ~ NA, .default = rachisLength),
         stemDiameterUpper = case_when(stemDiameterUpper < 5 | stemDiameterUpper > 16 ~ NA, .default = stemDiameterUpper), 
         stemDiameterLower = case_when(stemDiameterLower < 12.5 | stemDiameterLower > 33 ~ NA, .default = stemDiameterLower),
         tillersPerPlant = case_when(tillersPerPlant > 4.5 ~ NA, .default = tillersPerPlant), 
         thirdLeafWidth = case_when(thirdLeafWidth < 3 ~ NA, .default = thirdLeafWidth),
         thirdLeafLength = case_when(thirdLeafLength > 90 ~ NA, .default = thirdLeafLength),
         plantHeight = case_when(plantHeight > 200 ~ NA, .default = plantHeight), 
         flagLeafWidth = case_when(flagLeafWidth > 8.5 ~ NA, .default = flagLeafWidth), 
         flagLeafLength = case_when(flagLeafLength > 62 ~ NA, .default = flagLeafLength), 
         estimatedPlotYield = case_when(estimatedPlotYield > 975 ~ NA, .default = estimatedPlotYield), 
         panicleGrainWeight = case_when(panicleGrainWeight > 85 ~ NA, .default = panicleGrainWeight),
         paniclesPerPlot = case_when(paniclesPerPlot > 30 ~ NA, .default = paniclesPerPlot), 
         leafAngleStandardDeviation = case_when(leafAngleStandardDeviation > 15 ~ NA, .default = leafAngleStandardDeviation),
         medianLeafAngle = case_when(medianLeafAngle < 20 | medianLeafAngle > 75 ~ NA, .default = medianLeafAngle), 
         daysToFlower = case_when(daysToFlower < 51 | daysToFlower > 82 ~ NA, .default = daysToFlower))

# Fit environment-wise BLUEs
bluesByEnvironment <- tibble(genotype = NULL, year = NULL, treatment = NULL, BLUE = NULL, intercept = NULL, phenotype = NULL)
envDatasets <- list(sap2020LN, sap2020HN, sap2021LN, sap2021HN)

for (i in 1:length(envDatasets))
{
  df <- envDatasets[i][[1]]
  phenotypes <- phenotypeList[i][[1]]
  columnKnots <- floor(max(df$column, na.rm = TRUE)/2) + 1
  rowKnots <- floor(max(df$row, na.rm = TRUE)/2) + 1
  for (j in phenotypes)
  {
    if(is.null(df[[j]]) | sum(!is.na(df[[j]])) < 5)
    {
      next
    }
    print(paste0(envYears[i], ':', envNitrogen[i], ':', j))
    
    model <- SpATS(j, genotype = 'genotype', genotype.as.random = FALSE, spatial = ~ SAP(column, row, nseg = c(columnKnots, rowKnots)), data = df)
    plot.SpATS(model, main = paste0(envYears[i], ':', envNitrogen[i], ':', j))
    summary <- summary(model)$coeff
    intercept <- summary['Intercept']
    summary <- summary %>%
      as_tibble(rownames = 'genotype') %>% 
      filter(str_detect(genotype, 'PI')) %>%
      rename(BLUE = value) %>%
      mutate(year = envYears[i],
             treatment = envNitrogen[i], 
             intercept = intercept, 
             phenotype = j)
      bluesByEnvironment <- bind_rows(bluesByEnvironment, summary)
  }
}

bluesByEnvironment <- bluesByEnvironment %>%
  rowwise() %>%
  mutate(valueID = str_c(year, treatment, phenotype, sep = ':'))

nitrogenTreatments <- c('LN', 'HN')
phenotypes <- unique(bluesByEnvironment$phenotype)
years <- c('2020', '2021')

for(y in years)
{
  for(t in nitrogenTreatments)
  {
    for(p in phenotypes)
    {
      currID <- paste0(y, ':', t, ':', p)
      if(currID %in% unique(bluesByEnvironment$valueID))
      {
        df <- bluesByEnvironment %>%
          dplyr::filter(valueID==currID)
        printHistogram(df, 'BLUE', title = currID)
      }
    }
  }
}

# second stage extreme value removal
bluesByEnvironment <- bluesByEnvironment %>%
  rowwise() %>%
  mutate(BLUE = case_when((year=='2021' & treatment=='HN' & phenotype=='seedMassPerPlant' & (BLUE > 50)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='paniclesPerPlant' & (BLUE < -2.5| BLUE > -0.25)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='percentStarch' & (BLUE < -10 | BLUE > 5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='percentAsh' & (BLUE < -0.37 | BLUE > 0.75)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='percentOil' & (BLUE < -0.5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='percentProtein' & (BLUE < -7.5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='percentMoisture' & (BLUE < -1.25 | BLUE > 2.75)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='stemDiameterUpper' & (BLUE > 8)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='stemDiameterLower' & (BLUE > 14.5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='tillersPerPlant' & (BLUE > 0.5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='thirdLeafWidth' & (BLUE < -2.5)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='thirdLeafLength' & (BLUE < -10 | BLUE > 30)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='plantHeight' & (BLUE > 75)) ~ NA,
                          (year=='2021' & treatment=='HN' & phenotype=='daysToFlower' & (BLUE < -5 | BLUE > 30)) ~ NA,
                          (year=='2021' & treatment=='LN' & phenotype=='stemDiameterUpper' & (BLUE > 6.25)) ~ NA,
                          (year=='2021' & treatment=='LN' & phenotype=='stemDiameterLower' & (BLUE > 15)) ~ NA,
                          (year=='2021' & treatment=='LN' & phenotype=='tillersPerPlant' & (BLUE > 0)) ~ NA,
                          (year=='2021' & treatment=='LN' & phenotype=='plantHeight' & (BLUE > 50)) ~ NA,
                          (year=='2021' & treatment=='LN' & phenotype=='daysToFlower' & (BLUE < -2.5 | BLUE > 30)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='percentStarch' & (BLUE < -14)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='percentAsh' & (BLUE > 0.4)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='percentOil' & (BLUE < -2 | BLUE > 1)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='percentProtein' & (BLUE < -2.5 | BLUE > 6)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='percentMoisture' & (BLUE < -1.25 | BLUE > 2)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='branchInternodeLength' & (BLUE > 2.5)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='primaryBranchNumber' & (BLUE > 63)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='rachisDiameterUpper' & (BLUE > 2.5)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='stemDiameterLower' & (BLUE < -7.5)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='tillersPerPlant' & (BLUE > 3)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='thirdLeafWidth' & (BLUE < -1.25 | BLUE > 4.25)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='thirdLeafLength' & (BLUE < -15 | BLUE > 30)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='plantHeight' & (BLUE > 25)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='extantLeafNumber' & (BLUE > 7.5)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='flagLeafWidth' & (BLUE < -0.25 | BLUE > 5.25)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='flagLeafLength' & (BLUE < -10 | BLUE > 29)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='estimatedPlotYield' & (BLUE < -380)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='panicleGrainWeight' & (BLUE > 65)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='paniclesPerPlot' & (BLUE > 10)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='leafAngleStandardDeviation' & (BLUE > 10)) ~ NA,
                          (year=='2020' & treatment=='HN' & phenotype=='medianLeafAngle' & (BLUE > 25)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='branchInternodeLength' & (BLUE > 1.75)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='primaryBranchNumber' & (BLUE > 63)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='rachisDiameterUpper' & (BLUE > 2.75)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='rachisDiameterLower' & (BLUE > 6.25)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='stemDiameterUpper' & (BLUE > 7.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='stemDiameterLower' & (BLUE > 17.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='tillersPerPlant' & (BLUE > 2.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='thirdLeafWidth' & (BLUE < -1.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='thirdLeafLength' & (BLUE > 25)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='plantHeight' & (BLUE > 20)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='extantLeafNumber' & (BLUE < -2.5 | BLUE > 7.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='flagLeafWidth' & (BLUE  < -1 | BLUE > 4.25)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='flagLeafLength' & (BLUE < -7.5| BLUE > 30)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='estimatedPlotYield' & (BLUE > 250)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='leafAngleStandardDeviation' & (BLUE < -6 | BLUE > 4.5)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='medianLeafAngle' & (BLUE > 25)) ~ NA,
                          (year=='2020' & treatment=='LN' & phenotype=='daysToFlower' & (BLUE > 27.5)) ~ NA,
                          .default = BLUE))

# Calculate nitrogen response
bluesByEnvironmentWide <- bluesByEnvironment %>%
  select(!valueID) %>%
  rowwise() %>%
  mutate(expectedValue = BLUE + intercept) %>%
  pivot_wider(id_cols = c(genotype, year, phenotype), 
              names_from = treatment, values_from = expectedValue) %>%
  filter(!is.na(LN) & !is.na(HN)) %>%
  mutate(NR = (HN - LN)/LN, 
         valueID = str_c(year, phenotype, sep = ':'))

for(y in years)
{
    for(p in phenotypes)
    {
      currID <- paste0(y, ':', p)
      if(currID %in% unique(bluesByEnvironmentWide$valueID))
      {
        df <- bluesByEnvironmentWide %>%
          dplyr::filter(valueID==currID)
        printHistogram(df, 'NR', title = currID)
      }
    }
}

# Remove outliers
bluesByEnvironmentWide <- bluesByEnvironmentWide %>%
  mutate(NR = case_when((year=='2021' & phenotype=='stemDiameterUpper' & (NR < -0.25 | NR > 0.35)) ~ NA,
                        (year=='2021' & phenotype=='stemDiameterLower' & (NR > 0.15)) ~ NA,
                        (year=='2021' & phenotype=='tillersPerPlant' & (NR < -250)) ~ NA,
                        (year=='2021' & phenotype=='plantHeight' & (NR < -0.5 | NR > 0.5)) ~ NA,
                        (year=='2021' & phenotype=='daysToFlower' & (NR < -0.2 | NR > 0.15)) ~ NA,
                        (year=='2020' & phenotype=='branchInternodeLength' & (NR > 1.25)) ~ NA,
                        (year=='2020' & phenotype=='primaryBranchNumber' & (NR < -0.25 | NR > 0.625)) ~ NA,
                        (year=='2020' & phenotype=='rachisDiameterUpper' & (NR > 1.25)) ~ NA,
                        (year=='2020' & phenotype=='rachisDiameterLower' & (NR < -0.3 | NR > 0.6)) ~ NA,
                        (year=='2020' & phenotype=='rachisLength' & (NR < -0.4 | NR > 0.625)) ~ NA,
                        (year=='2020' & phenotype=='stemDiameterUpper' & (NR > 0.5)) ~ NA,
                        (year=='2020' & phenotype=='stemDiameterLower' & (NR < -0.3 | NR > 0.4)) ~ NA,
                        (year=='2020' & phenotype=='tillersPerPlant' & (NR  < -12.5 | NR > 12.5)) ~ NA,
                        (year=='2020' & phenotype=='thirdLeafWidth' & (NR < -0.375 | NR > 0.25)) ~ NA,
                        (year=='2020' & phenotype=='thirdLeafLength' & (NR > 0.4)) ~ NA,
                        (year=='2020' & phenotype=='plantHeight' & (NR > 0.5)) ~ NA,
                        (year=='2020' & phenotype=='extantLeafNumber' & (NR < -0.375 | NR > 0.3)) ~ NA,
                        (year=='2020' & phenotype=='flagLeafWidth' & (NR > 0.75)) ~ NA,
                        (year=='2020' & phenotype=='flagLeafLength' & (NR > 0.75)) ~ NA,
                        (year=='2020' & phenotype=='estimatedPlotYield' & (NR < -6 | NR > 6)) ~ NA,
                        (year=='2020' & phenotype=='panicleGrainWeight' & (NR > 5)) ~ NA,
                        (year=='2020' & phenotype=='paniclesPerPlot' & (NR > 2)) ~ NA,
                        (year=='2020' & phenotype=='leafAngleStandardDeviation' & (NR > 1.5)) ~ NA,
                        (year=='2020' & phenotype=='medianLeafAngle' & (NR > 0.55)) ~ NA,
                        (year=='2020' & phenotype=='daysToFlower' & (NR < -0.2 | NR > 0.1)) ~ NA,
                        .default = NR))

# Fit BLUEs across years
sap2020LN$year <- '2020'
sap2020HN$year <- '2020'
sap2021LN$year <- '2021'
sap2021HN$year <- '2021'
sapLN <- bind_rows(sap2020LN, sap2021LN) %>%
  mutate(year = as.factor(year),
         row = case_when(year=='2021' ~ row + 100, .default = row),
         column = case_when(year=='2021' ~ column + 100, .default = column))
sapHN <- bind_rows(sap2020HN, sap2021HN) %>%
  mutate(year = as.factor(year),
         row = case_when(year=='2021' ~ row + 100, .default = row),
         column = case_when(year=='2021' ~ column + 100, .default = column))
nitrogenDatasets <- list(sapLN, sapHN)
dfNitrogen <- c('LN', 'HN')
commonPhenotypes <- intersect(phenotypes2020, phenotypes2021)
commonPhenotypesLN <- c('daysToFlower', 'plantHeight', 'tillersPerPlant', 'stemDiameterUpper', 'stemDiameterLower')
phenotypeList <- list(commonPhenotypesLN, commonPhenotypes)
bluesAcrossYears <- tibble(genotype = NULL, treatment = NULL, BLUE = NULL, intercept = NULL, phenotype = NULL)


for (i in 1:length(nitrogenDatasets))
{
  df <- nitrogenDatasets[i][[1]]
  phenotypes <- phenotypeList[i][[1]]
  columnKnots <- 35
  rowKnots <- 35
  for (j in phenotypes)
  {
    if(is.null(df[[j]]) | sum(!is.na(df[[j]])) < 5)
    {
      next
    }
    print(paste0(dfNitrogen[i], ':', j))
    
    model <- SpATS(j, genotype = 'genotype', genotype.as.random = FALSE, fixed = ~ year, spatial = ~ SAP(column, row, nseg = c(columnKnots, rowKnots)), data = df)
    plot.SpATS(model, main = paste0(dfNitrogen[i], ':', j))
    summary <- summary(model)$coeff
    intercept <- summary['Intercept']
    summary <- summary %>%
      as_tibble(rownames = 'genotype') %>% 
      filter(str_detect(genotype, 'PI')) %>%
      rename(BLUE = value) %>%
      mutate(treatment = dfNitrogen[i], 
             intercept = intercept, 
             phenotype = j)
    bluesAcrossYears <- bind_rows(bluesAcrossYears, summary)
  }
}

bluesAcrossYears <- bluesAcrossYears %>%
  rowwise() %>%
  mutate(valueID = str_c(treatment, phenotype, sep = ':'))

for(t in nitrogenTreatments)
  {
    for(p in commonPhenotypes)
    {
      currID <- paste0(t, ':', p)
      if(currID %in% unique(bluesAcrossYears$valueID))
      {
        df <- bluesAcrossYears %>%
          dplyr::filter(valueID==currID)
        printHistogram(df, 'BLUE', title = currID)
      }
    }
}

# Second stage outlier removal
bluesAcrossYears <- bluesAcrossYears %>%
  mutate(BLUE = case_when((treatment=='HN' & phenotype=='percentStarch' & (BLUE < -12.5)) ~ NA, 
                          (treatment=='HN' & phenotype=='percentAsh' & (BLUE > 0.75)) ~ NA, 
                          (treatment=='HN' & phenotype=='percentOil' & (BLUE < 0)) ~ NA, 
                          (treatment=='HN' & phenotype=='percentProtein' & (BLUE < -6.25 | BLUE > 3.75)) ~ NA, 
                          (treatment=='HN' & phenotype=='percentMoisture' & (BLUE < -0.75 | BLUE > 2.5)) ~ NA, 
                          (treatment=='HN' & phenotype=='stemDiameterUpper' & (BLUE > 7.5)) ~ NA, 
                          (treatment=='HN' & phenotype=='tillersPerPlant' & (BLUE > 0.25)) ~ NA, 
                          (treatment=='HN' & phenotype=='thirdLeafWidth' & (BLUE < -1.25)) ~ NA, 
                          (treatment=='HN' & phenotype=='thirdLeafLength' & (BLUE > 30)) ~ NA, 
                          (treatment=='HN' & phenotype=='plantHeight' & (BLUE > 100)) ~ NA, 
                          (treatment=='HN' & phenotype=='extantLeafNumber' & (BLUE > 5.5)) ~ NA, 
                          (treatment=='HN' & phenotype=='daysToFlower' & (BLUE > 30)) ~ NA, 
                          (treatment=='LN' & phenotype=='tillersPerPlant' & (BLUE > 1)) ~ NA, 
                          (treatment=='LN' & phenotype=='plantHeight' & (BLUE > 40)) ~ NA, 
                          .default = BLUE))

# Calculate NR
bluesAcrossYearsWide <- bluesAcrossYears %>%
  select(!valueID) %>%
  rowwise() %>%
  mutate(expectedValue = BLUE + intercept) %>%
  pivot_wider(id_cols = c(genotype, phenotype), 
              names_from = treatment, values_from = expectedValue) %>%
  filter(!is.na(LN) & !is.na(HN)) %>%
  mutate(NR = (HN - LN)/LN, 
         valueID = str_c(phenotype, sep = ':'))

for(p in commonPhenotypes)
  {
    currID <- p
    if(currID %in% unique(bluesAcrossYearsWide$valueID))
    {
      df <- bluesAcrossYearsWide %>%
        dplyr::filter(valueID==currID)
      printHistogram(df, 'NR', title = currID)
    }
}

# Remove outliers
bluesAcrossYearsWide <- bluesAcrossYearsWide %>% 
  mutate(NR = case_when((phenotype=='stemDiameterUpper' & (NR > 0.3)) ~ NA, 
                        (phenotype=='stemDiameterLower' & (NR < -0.1 | NR > 0.25)) ~ NA, 
                        (phenotype=='tillersPerPlant' & (NR < -0.25 | NR > 1.5)) ~ NA, 
                        (phenotype=='plantHeight' & (NR < -0.4 | NR > 0.125)) ~ NA, 
                        (phenotype=='daysToFlower' & (NR < -0.2 | NR > 0.15)) ~ NA, 
                        .default = NR))

# Fit overall BLUEs
df <- bind_rows(sapLN, sapHN) %>%
  mutate(treatment = as.factor(treatment))
bluesOverall <- tibble(genotype = NULL, BLUE = NULL, intercept = NULL, phenotype = NULL)
columnKnots <- 35
rowKnots <- 35
for (j in commonPhenotypesLN)
{
  if(is.null(df[[j]]) | sum(!is.na(df[[j]])) < 5)
  {
    next
  }
  print(j)
  
  model <- SpATS(j, genotype = 'genotype', genotype.as.random = FALSE, fixed = ~ year + treatment, spatial = ~ SAP(column, row, nseg = c(columnKnots, rowKnots)), data = df)
  plot.SpATS(model, main = j)
  summary <- summary(model)$coeff
  intercept <- summary['Intercept']
  summary <- summary %>%
    as_tibble(rownames = 'genotype') %>% 
    filter(str_detect(genotype, 'PI')) %>%
    rename(BLUE = value) %>%
    mutate(intercept = intercept, 
           phenotype = j)
  bluesOverall <- bind_rows(bluesOverall, summary)
}

for(p in commonPhenotypesLN)
{
  tmp <- filter(bluesOverall, phenotype==p)
  printHistogram(tmp, 'BLUE', p)
}

bluesOverall <- bluesOverall %>%
  mutate(BLUE = case_when((phenotype=='stemDiameterLower' & (BLUE < -2.5)) ~ NA, 
                          (phenotype=='tillersPerPlant' & (BLUE > 0)) ~ NA, 
                          (phenotype=='plantHeight' & (BLUE > 75)) ~ NA, 
                          .default = BLUE))
# Aggregate all BLUEs
bluesByEnvironmentWide <- bluesByEnvironmentWide %>%
  rename(BLUE = NR) %>%
  mutate(valueID = str_c(valueID, 'NR', sep = ':')) %>%
  select(genotype, BLUE, valueID)

bluesAcrossYearsWide <- bluesAcrossYearsWide %>%
  rename(BLUE = NR) %>%
  mutate(valueID = str_c(valueID, 'NR', sep = ':'))

bluesOverall <- bluesOverall %>%
  mutate(valueID = phenotype)

blues <- bind_rows(bluesByEnvironment, bluesByEnvironmentWide, bluesAcrossYears, bluesAcrossYearsWide, bluesOverall) %>%
  select(genotype, BLUE, valueID) %>% 
  pivot_wider(id_cols = genotype, values_from = BLUE, names_from = valueID)

# export BLUES 
write.csv(blues, 'BLUEs_SAP2020_2021.csv', quote = FALSE, row.names = FALSE)

# unify years of data with and without extreme values
sap2020$year <- '2020'
sap2021$year <- '2021'
sapUnifiedRaw <- bind_rows(sap2020, sap2021)

write.csv(sapUnifiedRaw, 'SAP_2020_2021_RAW.csv', quote = FALSE, row.names = FALSE)
write.csv(df, 'SAP_2020_2021_LowConfidenceAndExtremeValuesRemoved.csv', quote = FALSE, row.names = FALSE)
