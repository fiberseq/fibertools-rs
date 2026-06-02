# `fibertools-rs` end-to-end benchmarks

Criterion-driven wall-clock benchmarks for the four read-heavy fibertools
subcommands. Each iterates or emits per-record annotations (m6a / cpg /
msp / nuc / fire), so wall-clock is dominated by annotation traversal
cost rather than ML inference or external I/O:

- `pileup` (with `--m6a --cpg` to exercise every annotation type)
- `fire` (the only bench that writes annotations back to an output BAM)
- `track-decorators` (read-only traversal building a BED12 decorator track)
- `extract --all` (per-record TSV with every annotation column)

Every bench uses `sample_size(10)` and otherwise leaves criterion's timing
defaults alone. Slower benches just take longer to gather their 10 samples;
criterion may print an "unable to complete 10 samples" note, which is
harmless. A full `cargo bench --bench subcommands` run takes ~15–20 minutes
on a warm machine.

## 1. Fetch the dataset (once)

```sh
bash benches/fetch-data.sh
```

Pulls ~133 MB of public FIRE test-data into `benches/data/fire-test-data/`
(gitignored). Idempotent — re-running skips files already present. No
AWS CLI required; the bucket is anonymous-public over plain HTTPS.

## 2. Save a `main` baseline

On the commit you want to compare against (typically `main` HEAD):

```sh
cargo bench --bench subcommands -- --save-baseline main
```

Persists per-bench measurements to `target/criterion/<group>/main/`.
`target/criterion/` lives outside git, so the baseline survives any
`git switch` / `git rebase`.

## 3. Compare a feature branch

After switching to the branch you want to measure:

```sh
cargo bench --bench subcommands -- --baseline main
```

Criterion prints a per-bench table:

```
pileup/default          time:   [12.341 s 12.567 s 12.812 s]
                        change: [+3.2% +5.1% +7.0%] (p = 0.001 < 0.05)
                        Performance has regressed.
```

HTML reports land in `target/criterion/<group>/report/`.

## Running individual benches

```sh
cargo bench --bench subcommands -- pileup        # only the pileup group
cargo bench --bench subcommands -- --quick fire  # one sample, fast smoke
cargo bench --bench subcommands -- extract       # extract group only
```

## Reducing measurement noise

Wall-clock benches are sensitive to background load. Recommended:

- **Linux:** set the CPU governor to `performance` for the duration of the run
  (`sudo cpupower frequency-set -g performance`).
- **macOS:** plug into AC power, close Chrome and other heavyweights,
  consider excluding `benches/data/` from Spotlight (System Settings →
  Siri & Spotlight → Spotlight Privacy).
- Run with `nice -n -5 cargo bench ...` if you have the permission to
  raise priority, and avoid concurrent `cargo build` in another shell.

Variance >5% between back-to-back runs of the same commit means the
environment is too noisy to trust deltas under ~10%; address the noise
before reading the comparison.

## Refreshing the baseline

If `main` moves and you want a new comparison reference:

```sh
git switch main && git pull
cargo bench --bench subcommands -- --save-baseline main
git switch -                              # back to your branch
cargo bench --bench subcommands -- --baseline main
```

Criterion overwrites the `main` baseline in place; no cleanup needed.

## Why these four subcommands?

They each do substantial per-record annotation work and none of them
are dominated by orthogonal cost centers (ML inference, pangenome-graph
ops, network I/O) that would swamp the annotation-traversal signal:

- `pileup` — builds per-position coverage from m6a/cpg/msp/nuc/fire annotations
- `fire` — reads, scores, and writes annotations on every fiber
- `track-decorators` — emits a BED12 plus a decorator track from fire/m6a/cpg/nuc positions
- `extract --all` — dumps every per-record annotation column to TSV

`predict_m6a`, `add_nucleosomes`, and `pg_*` are intentionally not
benched here — their wall-clock is dominated by ML inference or
pangenome-graph operations, which would dilute the signal we care about.
