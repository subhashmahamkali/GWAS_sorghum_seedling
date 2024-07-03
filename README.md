
# Sorghum plasicity: trade-off between growth and reproduction under temperature and N stressors

Plant plasticity or resilience requires the efficient use of available resources to fuel growth while maintaining function. 
Environmental stressors, such as high temperatures and nitrogen limitation, can challenge the balance between growth and reproduction. 
In this study, we are leveraging existing datasets in sorghum to investigate phenotypic plasticity under temperature and nitrogen stresses.


## Genomics analysis of Sorghum Assocation Panel (SAP) 

#### Sorghum SNP Calling on V5:


### old version RAW variants:-43811787
compared plink2.0 by calculating these basic stats manually result was corelated.
now filtering these vcf with maf=0.05 and missing rate=0.7 using plink2.0



# Phenotype


# Quantgen analysis

### Seedling phenotypes

we have collected sorghum seedling data from UNL greenhouse grown under both high and low nitrogen conditions.
we have grown 346 genotupes from SAP(sorghum association panel with 1 check).
The data collected are plant height, leaf count, shoot dry weight and shoot fresh weight and stored under DATA folder 
The .RMD files in the profiling are used to calculate BLUPs and Phenotype distribution with calculated BLUP values.

Here entire main data and sub-data sets are in the data folder.
Blup_out files are saved under the output folder.
all the code is in the profiling folder.
after creating the histograms and density plots they are saved under the graphs folder.

update 03.20.2024

curently working on genotype data collected from WGS of SAP

using plink2.0 to estimate basic stats (minor allele freequecny and misssing rate) from vcf file (RAW) 


