# Benchmark results: `molecular-annotation` vs `main`

Criterion wall-clock means, comparing the `base` baseline (main, saved 2026-05-28)
against the `molecular-annotation` baseline (this PR, saved 2026-06-10).

| Benchmark | main (`base`) | this PR (`molecular-annotation`) | Change  |
|-----------|--------------:|---------------------------------:|--------:|
| decorator |       3.729 s |                          3.416 s | −8.4%   |
| extract   |       4.919 s |                          4.692 s | −4.6%   |
| fire      |       7.859 s |                          6.699 s | −14.8%  |
| pileup    |       8.104 s |                          7.789 s | −3.9%   |

Every subcommand is faster on this branch, with `fire` the biggest win (~15%).

**Note:** ignore the `−39%` that criterion's built-in `pileup/change/` folder reports —
it was computed against an anomalous `new` baseline (4.75 s, a different/interrupted
workload) rather than the clean `base` snapshot. The honest pileup number is −3.9%.

_Source: `target/criterion/<group>/{base,molecular-annotation}/estimates.json`._
