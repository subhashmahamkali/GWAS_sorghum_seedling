#!/bin/sh
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --job-name=picard_PI656082
#SBATCH --mail-user=smahamkalivenkatas2@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/work/jyanglab/subhash/variant_calling/4.picard/0.log_files/PI656082_picard.err
#SBATCH --output=/work/jyanglab/subhash/variant_calling/4.picard/0.log_files/PI656082_picard.out

ml picard/3.0
picard  MarkDuplicates \
INPUT=/work/jyanglab/subhash/variant_calling/3.alignment/PI656082.srt.bam \
OUTPUT=/work/jyanglab/subhash/variant_calling/4.picard/PI656082_picard_dedup.bam \
METRICS_FILE=/work/jyanglab/subhash/variant_calling/4.picard/PI656082_picard_metrics.txt \
CREATE_INDEX=true \
REMOVE_D
REMOVE_DUPLICATES=true

