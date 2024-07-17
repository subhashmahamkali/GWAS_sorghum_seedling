#### Sorghum SNP Calling on V5:

- Raw input data (`trimmed fastq`) were copied from Hongyu: `/work/jyanglab/subhash/variant_calling/1.trimmed_data/`
- Data folders: `/work/jyanglab/subhash/variant_calling`
- Code: `profiling/1.variant_calling/`

- Important files: 
  - reference genome: `/work/jyanglab/subhash/variant_calling/2.refrence/Phytozome/PhytozomeV13/Sbicolor/v5.1/assembly/Sbicolor_730_v5.0.fa`
  - Deduplicated BAM files: `/work/jyanglab/subhash/variant_calling/4.picard/` (n: 400 bam files)
  - gVCF file: `/work/jyanglab/subhash/variant_calling/5.gvcf/*.g.vcf`
  - Raw VCF file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/3.merged_vcf/RAW_SAP_BQSR.vcf.gz` (116G)
  - A filtered VCF file: `/work/jyanglab/subhash/variant_calling/9.gvcf_to_vcf/4.filtering/1.annotate/SAP_V5_annotate.vcf` (289G)
  
  
#### Overall pipeline for varinat calling 

##### Step 1: fastqc and trimming the raw reads

Here I have used fastp to trim the low quality reads and adapter sequences from RAW unmapped reads
file path - 

```
# quality check
fastqc -o /work/jyanglab/subhash/variant_calling/1.trimmed_data/0.1_fastqc_report/ /work/jyanglab/subhash/variant_calling/1.trimmed_data/{PI}_2.fastq.gz

# trim the low quality reads
ref="/work/jyanglab/subhash/variant_calling/2.refrence/Phytozome/PhytozomeV13/Sbicolor/v5.1/assembly/Sbicolor_730_v5.0.fa"
reads="/work/jyanglab/subhash/variant_calling/1.trimmed_data"
out="/work/jyanglab/subhash/variant_calling/3.alignment"
bwa mem -t 8 -R "@RG\\tID:{PI}\\tPL:ILLUMINA\\tSM:{PI}" ${ref} ${reads}/{PI}_1.fastq.gz ${reads}/{PI}_2.fastq.gz|samtools view --threads 8 -bS -f 2 -q 30 -h - | samtools sort --threads 8 -o ${out}/{PI}.srt.bam -'
```


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