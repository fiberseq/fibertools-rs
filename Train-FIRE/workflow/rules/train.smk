rule train_model:
    """Train one FIRE model on its region_set's training-data.bed.gz."""
    input:
        training=lambda wc: f"results/region_sets/{exp_region_set(wc.exp)}/training-data.bed.gz",
    output:
        gbdt="results/experiments/{exp}/FIRE.gbdt.json",
        xgb="results/experiments/{exp}/FIRE.xgb.json",
        conf="results/experiments/{exp}/FIRE.conf.json",
        fdr_pdf="results/experiments/{exp}/FIRE.FDR.pdf",
        feat_pdf="results/experiments/{exp}/FIRE.feature.importance.pdf",
        metrics="results/experiments/{exp}/metrics.json",
    conda:
        "../envs/env.yml"
    threads: 96
    resources:
        mem_mb=get_mem_mb,
    params:
        outdir=lambda wc: f"results/experiments/{wc.exp}",
        p=lambda wc: exp_train_params(wc.exp),
        script=workflow.source_path("../scripts/train-fire-model.py"),
        outer_jobs=12,
        inner_jobs=8,
    shell:
        r"""
        export OMP_NUM_THREADS={params.inner_jobs}
        python {params.script} \
          {input.training} \
          --outdir {params.outdir} \
          --outer-jobs {params.outer_jobs} \
          --inner-jobs {params.inner_jobs} \
          --train-fdr {params.p[train_fdr]} \
          --test-fdr {params.p[test_fdr]} \
          --subset-max-train {params.p[subset_max_train]} \
          --direction {params.p[direction]} \
          --min-msp-length-for-positive-fire-call {params.p[min_msp_length_for_positive_fire_call]} \
          --min-msp-length-for-negative-fire-call {params.p[min_msp_length_for_negative_fire_call]} \
          $( [ "{params.p[grid_search]}" = "True" ] && echo --grid-search ) \
          --n-estimators-grid "{params.p[n_estimators]}" \
          --max-depth-grid "{params.p[max_depth]}" \
          --min-child-weight-fracs "{params.p[min_child_weight_fracs]}" \
          --colsample-bytree-grid "{params.p[colsample_bytree]}" \
          --gamma-grid "{params.p[gamma]}" \
          --learning-rate-grid "{params.p[learning_rate]}" \
          --early-stopping-rounds {params.p[early_stopping_rounds]} \
          --early-stopping-val-frac {params.p[early_stopping_val_frac]} \
          --mokapot-max-iter {params.p[mokapot_max_iter]} \
          $( [ "{params.p[balance_train]}" = "True" ] && echo --balance-train || echo --no-balance-train ) \
          $( [ "{params.p[mokapot_override]}" = "True" ] && echo --mokapot-override )
        """


rule compare_models:
    """Aggregate per-experiment FDR curves + metrics into a single plot + TSV."""
    input:
        metrics=expand("results/experiments/{exp}/metrics.json", exp=EXPERIMENTS),
        confs=expand("results/experiments/{exp}/FIRE.conf.json", exp=EXPERIMENTS),
    output:
        pdf="results/comparison/fdr_overlay.pdf",
        tsv="results/comparison/metrics.tsv",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=get_mem_mb,
    params:
        exps=EXPERIMENTS,
        script=workflow.source_path("../scripts/compare_models.py"),
    shell:
        r"""
        python {params.script} \
          --root results/experiments \
          --experiments {params.exps} \
          --out-pdf {output.pdf} \
          --out-tsv {output.tsv}
        """
