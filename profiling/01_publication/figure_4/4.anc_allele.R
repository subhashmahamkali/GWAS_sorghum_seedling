# Check ancestral allele at the lead SNP

awk '$1==9 && $2==59161829' /work/jyanglab/subhash/sorgsd/bal_s/ans_freq/Sorghum_ancestral_allele_V3.1.txt

#9       59161829        G       0       A       1       A

# Get allele frequencies for Wild
bcftools view \
-S /work/jyanglab/subhash/sorgsd/xpclr/sample_lists/wild.clean.txt \
-r 9:59156829-59166829 \
/work/jyanglab/subhash/sorgsd/filterfed_vcf/SorGSD.289snp.miss05.vcf.gz | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' | \
awk '$2==59160864 || $2==59161888 || $2==59161987 || $2==59162172 || $2==59162364 || $2==59163689 || $2==59163950'


# Improved
bcftools view \
-S /work/jyanglab/subhash/sorgsd/xpclr/sample_lists/improved.clean.txt \
-r 9:59156829-59166829 \
/work/jyanglab/subhash/sorgsd/filterfed_vcf/SorGSD.289snp.miss05.vcf.gz | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n'

# Landrace
bcftools view \
-S /work/jyanglab/subhash/sorgsd/xpclr/sample_lists/landrace.clean.txt \
-r 9:59156829-59166829 \
/work/jyanglab/subhash/sorgsd/filterfed_vcf/SorGSD.289snp.miss05.vcf.gz | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n'