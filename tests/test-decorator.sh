#!/bin/bash
set -euo pipefail

#scp hyak:'projects/phased-fdr-and-peaks/Train-FIRE-v2.0/FIRE.*.json' models/.

#samtools view -F 2308 -@ 16 -b chr20_PS00388_COLO829BL_1.bam chr20:0-1000000 |
samtools view -F 2308 -@ 16 -b tests/data/chr20.hifi.bam chr20:0-1000000 |
    cargo run --release -- fire - xx.tmp.bam

#
cargo run --release -- track xx.tmp.bam --bed12 tmp.bed12.gz | sort -k1,1 -k2,2n | bgzip -@ 8 >tmp.dec.bed12.gz

rm xx.tmp.bam

# samtools view -F 16 -@ 16 -b - |

echo ""
rg FIRE -z tmp.dec.bed12.gz | cut -f 9,10 | sort | datamash -g 1 sum 2
echo ""

as="../FIRE/workflow/templates/bed12_filter.as"
bedToBigBed -allow1bpOverlap -type=bed12+ -as=$as tmp.bed12.gz tests/data/hg38.analysisSet.fa.fai tmp.bed12.bb

echo "done bed12"

as="../FIRE/workflow/templates/decoration.as"
bedToBigBed -allow1bpOverlap -type=bed12+ -as=$as tmp.dec.bed12.gz tests/data/hg38.analysisSet.fa.fai tmp.dec.bed12.bb

aws s3 sync --profile Mitchell_Vollger --exclude "*" --include "tmp.*bb" $@ . s3://stergachis-public1/Mitchell/temp/FIREv2/

aws s3 sync --profile Mitchell_Vollger --exclude "*" --include "hub.txt" $@ . s3://stergachis-public1/Mitchell/temp/FIREv2/
