rule build_trackhub:
    input:
        base=expand("results/experiments/{exp}/fire-fibers.bb", exp=EXPERIMENTS),
        dec=expand(
            "results/experiments/{exp}/fire-fiber-decorators.bb", exp=EXPERIMENTS
        ),
        sizes="results/shared/chrom.sizes",
    output:
        hub="results/trackhub/hub.txt",
        genomes="results/trackhub/genomes.txt",
        trackdb="results/trackhub/trackDb.txt",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        name=config["trackhub"]["name"],
        short=config["trackhub"]["short_label"],
        long=config["trackhub"]["long_label"],
        email=config["trackhub"]["email"],
        genome=GENOME,
        exps=EXPERIMENTS,
    shell:
        r"""
        python workflow/scripts/make_trackhub.py \
          --hub-dir results/trackhub \
          --name "{params.name}" \
          --short-label "{params.short}" \
          --long-label "{params.long}" \
          --email "{params.email}" \
          --genome {params.genome} \
          --experiments {params.exps} \
          --results-root results/experiments \
          --chrom-sizes {input.sizes}
        """
