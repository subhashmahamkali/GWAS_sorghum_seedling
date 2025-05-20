#!/bin/sh
#SBATCH --job-name=filter_missing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=10G
#SBATCH --error=filter_missing.err
#SBATCH --output=filter_missing.out
#SBATCH --licenses=common
vcftools --vcf SorGSD.289snp.imputed.vcf --keep propinqaum_sample.txt --recode --out propinqaum_3sample
vcftools --vcf propinqaum_3sample.recode.vcf --freq --out propinqaum_3sample
vcftools --vcf SorGSD.289snp.imputed.vcf --keep landrace_samples.txt --recode --out landrace_107samples
