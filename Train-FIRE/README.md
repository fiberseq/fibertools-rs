# Train-FIRE — model training sweep

Snakemake workflow (following [SmkTemplate](https://github.com/mrvollger/SmkTemplate))
for training multiple FIRE models under different training configurations and
visualizing them together in a single UCSC track hub.

## Layout

```
config/config.yaml                shared inputs + region_sets + experiment list
workflow/Snakefile                entrypoint
workflow/rules/*.smk              region building, features, training, hub
workflow/scripts/*.py             training + aggregation + hub assembly
workflow/envs/env.yml             conda env for every rule
workflow/profiles/default/*.yaml  snakemake profile
resources/mixed-positives/        input peak BEDs
results/shared/                   shared sampled BAM + feature table
results/region_sets/<rs>/         positives/negatives/training-data (shared across experiments)
results/experiments/<exp>/        per-experiment trained models + FIRE calls
results/trackhub/                 UCSC track hub (hub.txt / trackDb.txt / bb/)
results/comparison/               FDR overlay plot + metrics TSV
```

## Quickstart

```bash
# from Train-FIRE/
pixi install
pixi run test                         # dry run
pixi run snakemake -j 32 --use-conda  # full run
pixi run snakemake -j 32 train_only   # train+compare, skip trackhub
```

## Adding an experiment

Edit [config/config.yaml](config/config.yaml) under `experiments:`. Each
experiment picks its own positives, negative strategy (`shuffle` or
`complement`), and optional `train:` overrides for the XGBoost grid search.

```yaml
experiments:
  my_new_model:
    positive_beds:
      - path: resources/mixed-positives/ATAC.bed.gz
        awk_filter: "$5 >= 1500"
      - path: resources/mixed-positives/peaks_CTCF_ENCFF951PEM.bed.gz
    negative_strategy: shuffle
    train:
      n_estimators: [300]
      max_depth: [12]
```

Feature extraction is shared: `ft fire -f` runs once over the union of every
experiment's positive + negative regions, and each experiment labels that
shared table via `bedtools intersect`.

## Test region

`test_region` in the config (default `chr19`) is held out from training
(injected into the exclusion bed automatically) and used as the sole region
for applying every trained model + trackhub visualization.

## Track hub

`results/trackhub/` is a standalone UCSC Track Hub. `rsync` it to a
web-accessible location and point the UCSC Genome Browser at `<url>/hub.txt`.
One bigBed9 FIRE-element track per experiment, visibility `squish`, colored
per model.
