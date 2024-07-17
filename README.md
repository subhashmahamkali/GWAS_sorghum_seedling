[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Sorghum plasicity: trade-off between growth and reproduction under temperature and N stressors

Plant plasticity or resilience requires the efficient use of available resources to fuel growth while maintaining function. 
Environmental stressors, such as high temperatures and nitrogen limitation, can challenge the balance between growth and reproduction. 
In this study, we are leveraging existing datasets in sorghum to investigate phenotypic plasticity under temperature and nitrogen stresses.


## Genomics and phenomics analyses of Sorghum Association Panel (SAP) 

See SNP calling pipeline [here](pipeline4snp.md)


#### Phenotype

##### Seedling phenotypes

we have collected sorghum seedling data from UNL greenhouse grown under both high and low nitrogen conditions.
we have grown 346 genotupes from SAP(sorghum association panel with 1 check).
The data collected are plant height, leaf count, shoot dry weight and shoot fresh weight and stored under DATA folder 
The .RMD files in the profiling are used to calculate BLUPs and Phenotype distribution with calculated BLUP values.

Here entire main data and sub-data sets are in the data folder.
Blup_out files are saved under the output folder.
all the code is in the profiling folder.
after creating the histograms and density plots they are saved under the graphs folder.


#### Quantgen analysis



### old version RAW variants:-43811787
compared plink2.0 by calculating these basic stats manually result was corelated.
now filtering these vcf with maf=0.05 and missing rate=0.7 using plink2.0

update 03.20.2024

curently working on genotype data collected from WGS of SAP

using plink2.0 to estimate basic stats (minor allele freequecny and misssing rate) from vcf file (RAW) 


------------

# Project Guideline

- To guide group members having a better sense about the project layout, here we briefly introduce the specific purposes of the [dir system](https://jyanglab.github.io/2017-01-07-project/). The layout of dirs is based on the idea borrowed from [ProjectTemplate](http://projecttemplate.net/architecture.html).

- The guideline for the collaborative [workflow](https://jyanglab.github.io/2017-01-10-project-using-github/).

- Throw ideas and check out progress and things [to-do](TODO.md). 

## License
This is an ongoing research project. It was intended for internal lab usage. It has not been extensively tested. Use at your own risk.
It is a free and open source software, licensed under [GPLv3](LICENSE).
