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

Experiments reference a named **region_set** (positives + negative-exclude
beds) defined under `region_sets:`. Every experiment sharing a region_set
shares the regions/features/training-data pipeline — only the XGBoost
training step runs per-experiment.

```yaml
region_sets:
  my_regions:
    positive_beds:
      - path: resources/mixed-positives/ATAC.bed.gz
        awk_filter: "$5 >= 1500"
      - path: resources/mixed-positives/peaks_CTCF_ENCFF951PEM.bed.gz
    negative_exclude_beds:
      - path: resources/mixed-positives/GM12878-fire-v0.1-peaks.bed.gz

experiments:
  my_new_model:
    region_set: my_regions
    train:
      n_estimators: [300]
      max_depth: [12]
```

Feature extraction is shared: `ft fire -f` runs once over the union of every
region_set's positive + negative regions, and each region_set labels that
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
