[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Sorghum plasicity: trade-off between growth and reproduction under temperature and N stressors

Plant plasticity or resilience requires the efficient use of available resources to fuel growth while maintaining function. 
Environmental stressors, such as high temperatures and nitrogen limitation, can challenge the balance between growth and reproduction. 
In this study, we are leveraging existing datasets in sorghum to investigate phenotypic plasticity under temperature and nitrogen stresses.


## Genomics and phenomics analyses of Sorghum Assocation Panel (SAP) 

#### Sorghum SNP Calling on V5:

- Data: HCC:`/work/jyanglab/subhash/variant_calling`
- Code: `1.variant_calling/`

- important files: 
  - reference genome: `/work/jyanglab/subhash/variant_calling/2.refrence/Phytozome/PhytozomeV13/Sbicolor/v5.1/assembly/Sbicolor_730_v5.0.fa`
  - Deduplicated BAM files: `/work/jyanglab/subhash/variant_calling/4.picard/` (n: 400 bam files)
  - Raw VCF file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/3.merged_vcf/RAW_SAP_BQSR.vcf.gz` (116G)
  - A filtered VCF file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/4.filtering/1.annotate/SAP_V5_annotate.vcf` (289G)
  
  
#### Overall pipeline for varinat calling 

##### Step 1: trimming the raw reads

Here I have used fastp to trim the low quality reads and adapter sequences from RAW unmapped reads
file path - 


##### Step 2: Aligning to reference genome

I have used BWA tool to align the trimmed reads to the reference genome ( version: 5)
then I used sam tools to convert into .bam format , then included only properly aloigned reads with mapping quality score of 30 or above and then sorted the aligned .bam file according to chromosome number and position.
path to reference genome-version-5: `/work/jyanglab/subhash/variant_calling/2.refrence/Phytozome/PhytozomeV13/Sbicolor/v5.1/assembly/Sbicolor_730_v5.0.fa`
path to aligned .bam files: `/work/jyanglab/subhash/variant_calling/3.alignment/*.srt.bam`


##### Step 3: de-duplication 

I have used picard tool to remove the duplicates ( which may arise during PCR amplification in library preparation). it also creates the index file for the .bam (for efficient quering and visualization)
path to the dedup .bam files: `/work/jyanglab/subhash/variant_calling/4.picard/*.bam`


##### Step 4: variant calling using GATK

here I have used haplotype caller from GATK (Genome Analysis Tool Kit) to perform variant calling on dedup .bam file, this outputs the g.vcf file ( this is genomic vcf file, which includes not only variants but also information about regions that are confidently non-variant)
This improves the accuracy and compelte coverage of variant sites and non variant sites.
path to gvcf: `/work/jyanglab/subhash/variant_calling/5.gvcf/*.g.vcf`


##### Step 5: joint genotyping  

- insted of calling variants for each sample which is more time and memory consuming process, I have called variants for every 5Mb interval across all the samples.

   5a) creating a data base for every 5Mb interval
   I have used GenomicsDBImport tool from GATK to import multiple GVCF files into 5Mb genomicsDB workspace for joint genotyping.
   path to DB:`/work/jyanglab/subhash/variant_calling/6.converting_to_vcf/1.database_5mb/`

   5b) converting g.vcf stored in database to vcf for every 5Mb interval
   I have used GenotypeGVCF tool from GATK to convert genomic variant data stored in GenomicsDB to VCF file. 
   path to vcf: `/work/jyanglab/subhash/variant_calling/6.converting_to_vcf/2.vcf_files/*.vcf.gz`

   5c) Merging all 5Mb vcf files into one vcf file.
   Here the input is multiple vcf files from 5Mb interval for the entire genome. I used GatherVCFs tool from GATK to merge into single VCF file.
   path to the RAW VCF ( also used as known sites for the downstream analysis: `/work/jyanglab/subhash/variant_calling/6.converting_to_vcf/3.merged_vcf/RAW_SAP.vcf.gz`


##### Step 6:  Base Quality Score Recalibration (BQSR)
I have used BaseRecalibrator tool from GATK to recalibrate each base quality score and stored this in a data table using the VCF file generated from above step as known sites. 
Then used Apply BQSR from GATK to apply the recalibration on dedu.bam files. Then, called the variants by follwing Step-4 and step-5 again and generated the final .vcf file.
path to BQSR g.vcf files: `/work/jyanglab/subhash/variant_calling/8.BQSR/3.g.vcf/*.g.vcf`
Path to final RAW  .vcf file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/3.merged_vcf/RAW_SAP_BQSR.vcf.gz`

##### Step 7:  filtered the VCF file 
After generating the RAW vcf file, I have used following parameters to filtering

```
mapping qaulity "MQ<20" 
Quality by Depth "QD<2.0"
Fisher Strand Bias "FS>60.0"
Mapping Quality Rank Sum Test "MQRankSum<-12.5"
Read Position Rank Sum Test "ReadPosRankSum<-8"
Path to filtered vcf file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/4.filtering/1.annotate/SAP_V5_annotate.vcf` 
```

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
