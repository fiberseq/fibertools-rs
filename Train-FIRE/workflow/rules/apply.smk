rule test_region_bam:
    """Extract test_region (default chr19) from the test alignment. Streams directly from s3 when given an s3:// URL."""
    input:
        aln=aln_input(TEST_BAM),
        fasta=FASTA,
    output:
        bam="results/shared/test_region.bam",
        csi="results/shared/test_region.bam.csi",
    conda:
        "../envs/env.yml"
    threads: 8
    resources:
        mem_mb=get_mem_mb,
    params:
        src=TEST_BAM,
        region=TEST_REGION,
    shell:
        r"""
        samtools view -@ {threads} -b --reference {input.fasta} \
            {params.src} {params.region} \
            -o {output.bam} --write-index
        """


rule apply_fire_bam:
    """Annotate the test-region BAM with FIRE calls from this experiment's model."""
    input:
        bam="results/shared/test_region.bam",
        model="results/experiments/{exp}/FIRE.gbdt.json",
        conf="results/experiments/{exp}/FIRE.conf.json",
    output:
        bam="results/experiments/{exp}/fire.bam",
        csi="results/experiments/{exp}/fire.bam.csi",
    conda:
        "../envs/env.yml"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        ft fire -t {threads} \
          --model {input.model} \
          --fdr-table {input.conf} \
          {input.bam} \
          | samtools sort -@ {threads} --write-index -o {output.bam}##idx##{output.csi} -
        """


def fire_bam_for_decorate(wildcards):
    """
    Source BAM for ft track-decorators:
      - baseline         -> the pre-annotated input CRAM subset to test_region
                            (FIRE calls from whatever model was baked into the input)
      - per-experiment   -> apply_fire_bam output with the freshly trained model
    """
    if wildcards.model == "baseline":
        return "results/shared/test_region.bam"
    return f"results/experiments/{wildcards.model}/fire.bam"


rule decorate_fibers:
    """Produce the base bed12 track + decorator overlay via `ft track-decorators`."""
    input:
        bam=fire_bam_for_decorate,
    output:
        base="results/models/{model}/fire-fibers.bed.gz",
        dec="results/models/{model}/fire-fiber-decorators.bed.gz",
        base_unsorted=temp("results/models/{model}/fire-fibers.unsorted.bed"),
    conda:
        "../envs/env.yml"
    threads: 8
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        samtools view -@ {threads} -u {input.bam} \
          | ft track-decorators -t {threads} --bed12 {output.base_unsorted} \
          | grep -v '^#' | grep -vw 'NUC' \
          | sort -k1,1 -k2,2n -k3,3n -k4,4 \
          | bgzip -@ {threads} > {output.dec}

        sort -k1,1 -k2,2n -k3,3n -k4,4 {output.base_unsorted} \
          | bgzip -@ {threads} > {output.base}
        """


rule base_to_bigbed:
    """bigBed 12+ for the FIRE-fibers base track."""
    input:
        bed="results/models/{model}/fire-fibers.bed.gz",
        sizes="results/shared/chrom.sizes",
        bed_as="workflow/templates/bed12_filter.as",
    output:
        bb="results/models/{model}/fire-fibers.bb",
        bed=temp("results/models/{model}/fire-fibers.bed"),
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        zcat -f {input.bed} > {output.bed}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={input.bed_as} \
          {output.bed} {input.sizes} {output.bb}
        """


rule decorator_to_bigbed:
    """bigBed 12+ decorator overlay track."""
    input:
        bed="results/models/{model}/fire-fiber-decorators.bed.gz",
        sizes="results/shared/chrom.sizes",
        dec_as="workflow/templates/decoration.as",
    output:
        bb="results/models/{model}/fire-fiber-decorators.bb",
        bed=temp("results/models/{model}/fire-fiber-decorators.bed"),
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        zcat -f {input.bed} > {output.bed}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={input.dec_as} \
          {output.bed} {input.sizes} {output.bb}
        """
