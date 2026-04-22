rule sample_bam:
    """Fractional subsample of the training alignment (local or remote URL) into a BAM."""
    input:
        aln=aln_input(TRAINING_BAM),
        fasta=FASTA,
    output:
        bam="results/shared/sample.bam",
        csi="results/shared/sample.bam.csi",
    conda:
        "../envs/env.yml"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
    params:
        src=TRAINING_BAM,
        rate=SAMPLE_RATE,
    shell:
        r"""
        samtools view -@ {threads} -b -s {params.rate} \
            --reference {input.fasta} \
            {params.src} \
            -o {output.bam} --write-index
        """


rule union_regions:
    """Concat every region_set's positives + negatives so we extract features from the BAM exactly once."""
    input:
        beds=all_rs_region_files(),
        fai=FAI,
    output:
        bed="results/shared/union_regions.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        zcat -f {input.beds} \
          | cut -f1-3 \
          | bedtools sort -g {input.fai} \
          | bedtools merge \
          | bgzip > {output.bed}
        """


rule regions_bam:
    """Subset the sampled BAM to reads overlapping any training region. Indexes
    the output so extract_features can scatter samtools views by chromosome."""
    input:
        bam="results/shared/sample.bam",
        regions="results/shared/union_regions.bed.gz",
    output:
        bam="results/shared/regions.bam",
        csi="results/shared/regions.bam.csi",
    conda:
        "../envs/env.yml"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        samtools view -@ {threads} -b -M \
            -L <(zcat -f {input.regions}) \
            {input.bam} -o {output.bam} --write-index
        """


rule extract_features:
    """Scatter ft fire -f by chromosome to work around the single BAM reader
    bottleneck inside a single ft process. Each chrom streams through its own
    samtools view | ft fire -f - pipeline; outputs are concatenated."""
    input:
        bam="results/shared/regions.bam",
        csi="results/shared/regions.bam.csi",
        regions="results/shared/union_regions.bed.gz",
    output:
        feats="results/shared/features.tsv.gz",
    conda:
        "../envs/env.yml"
    threads: 32
    resources:
        mem_mb=get_mem_mb,
    params:
        min_msp=TRAIN_DEFAULTS["min_msp_length_for_positive_fire_call"],
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        # Split the thread budget: outer = parallel chrom jobs, inner = ft threads.
        outer=$(( {threads} / 4 ))
        [ $outer -lt 1 ] && outer=1
        inner=4

        # Only scatter over chroms that actually appear in union_regions.
        zcat -f {input.regions} | cut -f1 | sort -u > "$tmpdir/chroms.txt"

        process_chrom() {{
            chrom="$1"
            samtools view -@ 1 -u "$BAM" "$chrom" \
              | ft fire -t "$INNER" --min-msp-length-for-positive-fire-call "$MIN_MSP" -f - \
              | awk 'NR == 1 || ($3 > $2 && $3 - $2 < 10000)' \
              > "$TMP/$chrom.tsv"
        }}
        export -f process_chrom
        export BAM={input.bam} MIN_MSP={params.min_msp} TMP="$tmpdir" INNER="$inner"

        xargs -a "$tmpdir/chroms.txt" -n1 -P "$outer" -I CHR bash -c 'process_chrom CHR'

        # Concat: take header from the first non-empty shard; rows from all.
        first=1
        for f in "$tmpdir"/*.tsv; do
            [ -s "$f" ] || continue
            if [ $first -eq 1 ]; then
                head -n 1 "$f"
                first=0
            fi
            tail -n +2 "$f"
        done | bgzip -@ 8 > {output.feats}
        """


rule build_training_data:
    """Label features: +1 if row overlaps positives (f>=0.25); -1 if overlaps negatives and no unfiltered positive."""
    input:
        feats="results/shared/features.tsv.gz",
        positives="results/region_sets/{rs}/positives.bed.gz",
        negatives="results/region_sets/{rs}/negatives.bed.gz",
        mask="results/region_sets/{rs}/neg_mask.bed.gz",
    output:
        bed="results/region_sets/{rs}/training-data.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        header=$(zcat -f {input.feats} | head -n 1 || true)
        (
          printf '%s\tLabel\n' "$header"
          bedtools intersect -f 0.25 -u -a {input.feats} -b {input.positives} \
            | sed 's/$/\t1/' | awk '/^chr/'
          bedtools intersect -f 0.25 -u -a {input.feats} -b {input.negatives} \
            | bedtools intersect -v -a - -b {input.mask} \
            | sed 's/$/\t-1/' | awk '/^chr/'
        ) | bgzip -@ 8 > {output.bed}
        echo "[{wildcards.rs}] training label counts:" >&2
        zcat -f {output.bed} | tail -n +2 | awk -F'\t' '{{print $NF}}' | sort | uniq -c >&2
        """
