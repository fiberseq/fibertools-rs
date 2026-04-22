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
    """Concat every experiment's positives + negatives so we extract features from the BAM exactly once."""
    input:
        beds=all_exp_region_files(),
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
    """Subset the sampled BAM to reads overlapping any training region."""
    input:
        bam="results/shared/sample.bam",
        regions="results/shared/union_regions.bed.gz",
    output:
        bam="results/shared/regions.bam",
    conda:
        "../envs/env.yml"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        samtools view -@ {threads} -b -M \
            -L <(zcat -f {input.regions}) \
            {input.bam} -o {output.bam}
        """


rule extract_features:
    """Run ft fire -f once to produce the shared MSP feature table."""
    input:
        bam="results/shared/regions.bam",
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
        ft fire -t {threads} --min-msp-length-for-positive-fire-call {params.min_msp} -f {input.bam} \
          | awk 'NR == 1 || ($3 > $2 && $3 - $2 < 10000)' \
          | bgzip -@ 8 > {output.feats}
        """


rule build_training_data:
    """Label features: +1 if row overlaps positives (f>=0.25); -1 if overlaps negatives and no unfiltered positive."""
    input:
        feats="results/shared/features.tsv.gz",
        positives="results/experiments/{exp}/positives.bed.gz",
        negatives="results/experiments/{exp}/negatives.bed.gz",
        mask="results/experiments/{exp}/neg_mask.bed.gz",
    output:
        bed="results/experiments/{exp}/training-data.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        (
          zcat -f {input.feats} | head -n 1 | sed 's/$/\tLabel/'
          bedtools intersect -f 0.25 -u -a {input.feats} -b {input.positives} \
            | sed 's/$/\t1/' | awk '/^chr/'
          bedtools intersect -f 0.25 -u -a {input.feats} -b {input.negatives} \
            | bedtools intersect -v -a - -b {input.mask} \
            | sed 's/$/\t-1/' | awk '/^chr/'
        ) | bgzip -@ 8 > {output.bed}
        echo "[{wildcards.exp}] training label counts:" >&2
        zcat -f {output.bed} | tail -n +2 | awk -F'\t' '{{print $NF}}' | sort | uniq -c >&2
        """
