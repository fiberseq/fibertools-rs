rule chrom_sizes:
    """UCSC chrom.sizes from the .fai (cols 1, 2). Shared across rules."""
    input:
        fai=FAI,
    output:
        sizes="results/shared/chrom.sizes",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"cut -f1,2 {input.fai} > {output.sizes}"


checkpoint build_chrom_list:
    """Snapshot the chroms that actually have reads in regions.bam so
extract_features can scatter per-chrom. Deriving from idxstats (vs the FAI)
avoids launching empty jobs for chroms with no training-region coverage."""
    input:
        bam="results/shared/regions.bam",
        csi="results/shared/regions.bam.csi",
    output:
        chroms="results/shared/chroms.txt",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        exclude_re=EXCLUDE_CHROMS_RE,
    shell:
        r"""
        samtools idxstats {input.bam} \
          | awk '$1 != "*" && $3 > 0 {{ print $1 }}' \
          | grep -Ev {params.exclude_re:q} > {output.chroms}
        """


rule fetch_reference:
    """Download hg38 analysis set from UCSC (once) and index it."""
    output:
        fa="resources/ref/hg38.analysisSet.fa",
        fai="resources/ref/hg38.analysisSet.fa.fai",
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
    shell:
        r"""
        curl -L --fail --retry 3 {params.url} \
          | gunzip -c > {output.fa}
        samtools faidx {output.fa}
        """
