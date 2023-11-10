#!/bin/bash
set -eu
bam=$1
mixed_positives=$2
negatives=$3

zcat -f $mixed_positives $negatives |
    shuf | head -n 20000 |
    sort -k1,1 -k2,2n | bedtools merge > \
    tmp.bed

samtools view -@ 16 -ML tmp.bed -b $bam -o tmp.bam --write-index

cargo run --release -- fire -t 16 -f tmp.bam | bgzip -@ 8 > \
    tmp.feats.gz

cat tmp.feats.gz | zcat | head -n 1 | sed 's/$/\tLabel/' > \
    traing-data.bed

bedtools intersect -u -a tmp.feats.gz -b $mixed_positives |
    sed 's/$/\t1/' >>traing-data.bed

bedtools intersect -u -a tmp.feats.gz -b $negatives |
    sed 's/$/\t-1/' >>traing-data.bed

bgzip -f -@ 16 traing-data.bed

rm tmp.bed tmp.bam tmp.bam.csi tmp.feats.gz
