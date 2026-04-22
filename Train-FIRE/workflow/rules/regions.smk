rule build_exclude:
    """Union of user-supplied exclusion beds + regex-excluded chroms + the held-out test chromosome."""
    input:
        beds=EXCLUDE_BEDS,
        fai=FAI,
    output:
        bed="results/shared/exclude.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        chrom_re=EXCLUDE_CHROMS_RE,
        test_region=TEST_REGION,
    shell:
        r"""
        ( \
            zcat -f {input.beds} | cut -f1-3; \
            awk -v OFS='\t' 'BEGIN{{re="{params.chrom_re}"}} $1 ~ re {{print $1,0,$2}}' {input.fai}; \
            awk -v OFS='\t' '$1 == "{params.test_region}" {{print $1,0,$2}}' {input.fai} \
        ) \
        | cut -f1-3 \
        | awk 'NR==FNR{{v[$1]; next}} ($1 in v)' {input.fai} - \
        | bedtools sort -g {input.fai} \
        | bedtools merge \
        | bgzip > {output.bed}
        """


rule build_positives:
    """Union region_set's positive beds (with optional awk_filter), merge, subsample to n_sites, subtract exclude."""
    input:
        beds=lambda wc: rs_positive_paths(wc.rs),
        exclude="results/shared/exclude.bed.gz",
        fai=FAI,
    output:
        bed="results/region_sets/{rs}/positives.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        n_sites=N_SITES,
        merge_dist=MERGE_DIST,
        sources=lambda wc: positive_source_cmds(wc.rs),
    shell:
        r"""
        ( {params.sources} ) \
          | awk 'NR==FNR{{v[$1]; next}} ($1 in v)' {input.fai} - \
          | bedtools sort -g {input.fai} \
          | bedtools merge -d {params.merge_dist} \
          | (shuf | head -n {params.n_sites} || true) \
          | bedtools sort -g {input.fai} \
          | bedtools subtract -a - -b {input.exclude} \
          | bgzip > {output.bed}
        printf "[{wildcards.rs}] positives: " >&2
        zcat -f {output.bed} | wc -l >&2
        """


rule build_neg_mask:
    """Union of positives + negative_exclude_beds; negatives will be kept off these regions."""
    input:
        beds=lambda wc: rs_neg_mask_paths(wc.rs),
        exclude="results/shared/exclude.bed.gz",
        fai=FAI,
    output:
        bed="results/region_sets/{rs}/neg_mask.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        merge_dist=MERGE_DIST,
    shell:
        r"""
        zcat -f {input.beds} \
          | cut -f1-3 \
          | awk 'NR==FNR{{v[$1]; next}} ($1 in v)' {input.fai} - \
          | bedtools sort -g {input.fai} \
          | bedtools merge -d {params.merge_dist} \
          | bedtools subtract -a - -b {input.exclude} \
          | bgzip > {output.bed}
        """


rule build_complement_negatives:
    """Everything not in the negative-exclusion mask, minus exclude."""
    input:
        mask="results/region_sets/{rs}/neg_mask.bed.gz",
        exclude="results/shared/exclude.bed.gz",
        fai=FAI,
    output:
        bed="results/region_sets/{rs}/complement_negatives.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        bedtools complement -i {input.mask} -g {input.fai} \
          | bedtools subtract -a - -b {input.exclude} \
          | bgzip > {output.bed}
        """


rule build_negatives:
    """Shuffle positives into the complement so negatives length-match positives on the same chromosome."""
    input:
        positives="results/region_sets/{rs}/positives.bed.gz",
        mask="results/region_sets/{rs}/neg_mask.bed.gz",
        complement="results/region_sets/{rs}/complement_negatives.bed.gz",
        fai=FAI,
    output:
        bed="results/region_sets/{rs}/negatives.bed.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        r"""
        bedtools shuffle \
            -excl {input.mask} \
            -incl {input.complement} \
            -i {input.positives} \
            -chrom -seed 42 -g {input.fai} \
          | sort -k1,1 -k2,2n \
          | bgzip > {output.bed}
        printf "[{wildcards.rs}] negatives: " >&2
        zcat -f {output.bed} | wc -l >&2
        """
