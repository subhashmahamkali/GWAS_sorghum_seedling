#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Sorghum V5 negative-ZSS genomic feature pipeline (with intergenic + thresholds)
# -----------------------------------------------------------------------------
# Run on HCC in directory containing:
#   zeroscore.bed
#   genomic_fetures/features.sorted.bed
#   genomic_fetures/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz
#
# Output:
#   genomic_fetures/feature_bp_exclusive_with_intergenic.tsv
#   zss_neg_feature_exclusive_intergenic_for_thresholds.tsv
# -----------------------------------------------------------------------------

BASE="/mnt/nrdstor/jyanglab/subhash/datasets/3.results/negative"
FEAT_DIR="${BASE}/genomic_fetures"
ZSS_BED="${BASE}/zeroscore.bed"
FEAT_BED="${FEAT_DIR}/features.sorted.bed"
GFF_GZ="${FEAT_DIR}/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz"

mkdir -p "${FEAT_DIR}/exclusive_with_intergenic"
WORK="${FEAT_DIR}/exclusive_with_intergenic"

echo "[1/8] Build negative ZSS BED (tab-delimited, sorted)..."
awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"}
     NF>=5 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $5<0 {
       print $1, int($2), int($3), int($4), $5, $1 ":" int($4)
     }' "${ZSS_BED}" \
| sort -k1,1 -k2,2n > "${BASE}/zss_neg.sorted.bed"

echo "[2/8] Build genome sizes from V5 GFF..."
gzip -dc "${GFF_GZ}" \
| awk 'BEGIN{OFS="\t"}
       /^##sequence-region/{
         c=$2;
         if(c ~ /^[0-9]+$/) c=sprintf("Chr%02d", c);
         if(c ~ /^Chr([0-9]|[0-9][0-9])$/) print c,$4
       }' \
| sort -k1,1V > "${WORK}/genome.sizes"

echo "[3/8] Merge intervals per feature class..."
for f in exon intron upstream_2kb downstream_2kb; do
  awk -v f="$f" 'BEGIN{FS=OFS="\t"} $4==f {print $1,$2,$3}' "${FEAT_BED}" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "${WORK}/${f}.merged.bed"
done

echo "[4/8] Create exclusive non-overlapping feature classes..."
cp "${WORK}/exon.merged.bed" "${WORK}/exon.excl.bed"

bedtools subtract -a "${WORK}/intron.merged.bed" \
                  -b "${WORK}/exon.excl.bed" \
> "${WORK}/intron.excl.bed"

cat "${WORK}/exon.excl.bed" "${WORK}/intron.excl.bed" \
| sort -k1,1 -k2,2n | bedtools merge -i - > "${WORK}/exon_intron.union.bed"

bedtools subtract -a "${WORK}/upstream_2kb.merged.bed" \
                  -b "${WORK}/exon_intron.union.bed" \
> "${WORK}/upstream_2kb.excl.bed"

cat "${WORK}/exon_intron.union.bed" "${WORK}/upstream_2kb.excl.bed" \
| sort -k1,1 -k2,2n | bedtools merge -i - > "${WORK}/exon_intron_upstream.union.bed"

bedtools subtract -a "${WORK}/downstream_2kb.merged.bed" \
                  -b "${WORK}/exon_intron_upstream.union.bed" \
> "${WORK}/downstream_2kb.excl.bed"

echo "[5/8] Build intergenic as complement of union(non-intergenic)..."
cat "${WORK}/exon.excl.bed" \
    "${WORK}/intron.excl.bed" \
    "${WORK}/upstream_2kb.excl.bed" \
    "${WORK}/downstream_2kb.excl.bed" \
| sort -k1,1 -k2,2n | bedtools merge -i - > "${WORK}/non_intergenic.union.bed"

bedtools complement -i "${WORK}/non_intergenic.union.bed" \
                    -g "${WORK}/genome.sizes" \
> "${WORK}/intergenic.excl.bed"

echo "[6/8] Build combined exclusive feature BED (with intergenic)..."
{
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"exon"}' "${WORK}/exon.excl.bed"
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"intron"}' "${WORK}/intron.excl.bed"
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"upstream_2kb"}' "${WORK}/upstream_2kb.excl.bed"
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"downstream_2kb"}' "${WORK}/downstream_2kb.excl.bed"
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"intergenic"}' "${WORK}/intergenic.excl.bed"
} | sort -k1,1 -k2,2n > "${WORK}/features_exclusive_with_intergenic.bed"

echo "[7/8] Intersect negative SNPs with exclusive feature set..."
bedtools intersect -a "${BASE}/zss_neg.sorted.bed" \
                   -b "${WORK}/features_exclusive_with_intergenic.bed" \
                   -wa -wb \
> "${BASE}/zss_neg_feature_exclusive_intergenic.tsv"

awk 'BEGIN{FS=OFS="\t"} {print $6,$5,$10}' \
  "${BASE}/zss_neg_feature_exclusive_intergenic.tsv" \
> "${BASE}/zss_neg_feature_exclusive_intergenic_for_thresholds.tsv"

echo "[8/8] Compute feature bp denominators..."
{
  awk 'BEGIN{s=0} {s+=($3-$2)} END{print "exon\t"s}' "${WORK}/exon.excl.bed"
  awk 'BEGIN{s=0} {s+=($3-$2)} END{print "intron\t"s}' "${WORK}/intron.excl.bed"
  awk 'BEGIN{s=0} {s+=($3-$2)} END{print "upstream_2kb\t"s}' "${WORK}/upstream_2kb.excl.bed"
  awk 'BEGIN{s=0} {s+=($3-$2)} END{print "downstream_2kb\t"s}' "${WORK}/downstream_2kb.excl.bed"
  awk 'BEGIN{s=0} {s+=($3-$2)} END{print "intergenic\t"s}' "${WORK}/intergenic.excl.bed"
} > "${FEAT_DIR}/feature_bp_exclusive_with_intergenic.tsv"

echo
echo "Done. Key outputs:"
echo "  ${BASE}/zss_neg_feature_exclusive_intergenic_for_thresholds.tsv"
echo "  ${FEAT_DIR}/feature_bp_exclusive_with_intergenic.tsv"
echo
wc -l "${BASE}/zss_neg.sorted.bed" "${BASE}/zss_neg_feature_exclusive_intergenic_for_thresholds.tsv"
cut -f3 "${BASE}/zss_neg_feature_exclusive_intergenic_for_thresholds.tsv" | sort | uniq -c
