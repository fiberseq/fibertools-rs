#!/usr/bin/env bash
set -euo pipefail

# Align PacBio HiFi reads to reference using minimap2
# Preserves all tags from input FASTA comments (-y)
# Uses extended CIGAR with =/X operators (--eqx)
# Outputs sorted CRAM

minimap2 -ax map-hifi -Y -y --eqx --MD ref.fa.gz reads.fa -p 1 \
    | samtools sort -O cram --reference ref.fa.gz -o test.cram -

samtools index test.cram
