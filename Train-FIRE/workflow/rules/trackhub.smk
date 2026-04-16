rule build_trackhub:
    input:
        base=expand("results/models/{model}/fire-fibers.bb", model=MODELS),
        dec=expand("results/models/{model}/fire-fiber-decorators.bb", model=MODELS),
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
        models=MODELS,
    shell:
        r"""
        python workflow/scripts/make_trackhub.py \
          --hub-dir results/trackhub \
          --name "{params.name}" \
          --short-label "{params.short}" \
          --long-label "{params.long}" \
          --email "{params.email}" \
          --genome {params.genome} \
          --models {params.models} \
          --results-root results/models \
          --chrom-sizes {input.sizes}
        """
