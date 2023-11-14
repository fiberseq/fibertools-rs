#!/bin/bash
set -eu
bam=$1
mixed_positives=./GM12878-mixed-positive-accessible.bed.gz
negatives=./GM12878-not-accessible.bed.gz

# zcat exclude.bed.gz GM12878-mixed-positive-accessible.bed.gz | rg -v _ | bedtools sort -g hg38.analysisSet.fa.fai | bedtools complement -g hg38.analysisSet.fa.fai -i - | rg -v _ | bgzip > GM12878-not-accessible.bed.gz

#zcat -f $mixed_positives $negatives |
#    shuf | head -n 20000 |
#    sort -k1,1 -k2,2n | bedtools merge > \
#    tmp.bed

zcat -f $mixed_positives | shuf | head -n 100000 | sort -k1,1 -k2,2n >tmp.pos.bed
bedtools shuffle -incl $negatives -i tmp.pos.bed -chrom -seed 42 -g hg38.analysisSet.fa.fai |
    sort -k1,1 -k2,2n >tmp.neg.bed
cat tmp.pos.bed tmp.neg.bed | sort -k1,1 -k2,2n >tmp.bed

echo "bed files done"

samtools view -@ 16 -b -ML tmp.bed $bam -o tmp.bam
cargo run --release -- fire -t 16 -f tmp.bam | sort -k1,1 -k2,2n | bgzip -@ 8 >tmp.feats.gz

echo "features done"

cat tmp.feats.gz | zcat | head -n 1 | sed 's/$/\tLabel/' >traing-data.bed

bedtools intersect -u -a tmp.feats.gz -b $mixed_positives | sed 's/$/\t1/' >>traing-data.bed

bedtools intersect -u -a tmp.feats.gz -b $negatives | sed 's/$/\t-1/' >>traing-data.bed

bgzip -f -@ 16 traing-data.bed

#rm tmp.*
