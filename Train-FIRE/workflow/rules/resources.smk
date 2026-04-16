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
