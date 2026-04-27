# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### 🚀 New subcommands

- **`ft call-peaks`**: peak caller for FIRE pileup BED output. Identifies enriched regions, applies optional FIRE fiber-level filters, and emits per-peak statistics.
- **`ft mock-fire`**: generate mock FIRE-style data for testing and benchmarking pipelines.
- **`ft benchmark`**: lightweight benchmarking harness for measuring throughput on key subcommands.

### 🔧 Improvements

- **`ft pileup`**: accepts multiple regions in a single invocation; default `--frac-fibers` filter applied to drop low-coverage windows.
- **FIRE fiber-level filters hoisted to global `FiberFilters`**: `--fire-filter`, `--min-msp`, `--min-ave-msp-size`, and `--skip-no-m6a` are now available as global flags shared by `ft fire`, `ft pileup`, and `ft call-peaks`. `--fire-filter` is a convenience preset; individual flags override its defaults.
- **Train-FIRE Snakemake pipeline** added under `Train-FIRE/` for end-to-end FIRE model training.
- Build/dev-loop improvements: thin-LTO and a dedicated test profile for faster iteration.

### 🐛 Bug fixes

- **fa:Z extras pairing on reverse-strand reads with overlapping/shared-start peaks (follow-up to #103)**: when multiple BED annotations share a query start position on a reverse-strand read, the post-hoc `ann_vals.reverse()` step in `FiberAnnotations::from_bam_tags` did not produce a permutation equivalent to the stable sort over flipped coordinates, scrambling the fa:Z extras (e.g. peak names/lengths) relative to their fs/fl annotations on extraction. Fix pre-pairs extras with annotations before the flip+sort so the permutation is consistent. BAM encodings were always correct; the bug only affected extraction. Affects v0.8.0–v0.8.2.

### 🧪 Tests

- `tests/call_peaks_test.rs`: end-to-end coverage of the new peak caller.
- `tests/fibertig_test.rs`: regression test for the fa:Z reverse-strand shared-start pairing bug, using Anna Minkina's exact failing input.
- Unit tests in `src/utils/input_bam.rs` covering `passes_fire_filter` semantics: skip-no-m6a behavior, min-msp, min-ave-msp-size, the `--fire-filter` combo, and explicit-flag overrides.

### ⚡ Performance

- Per-record FIRE feature extraction in `ft fire -f` is now parallelized (#104).

### [0.8.2]

#### Bug Fixes

- **`ft center`: fix off-by-one in `centered_query_start`/`centered_query_end` for minus strand reads.** The half-open interval negation was incorrect, causing m6a positions to index into the wrong base when slicing `query_sequence[-centered_query_end:]` on minus strand reads. Plus strand and `subset_sequence` were unaffected.
- **`ft center`: add missing `fire_qual` column to wide format header.** The header had 21 columns but data rows had 22, causing column misalignment when parsing.

### [0.8.1]

Fixed a bug where pg-inject could make fibertigs that were longer than needed at the end of chromosomes.

## [0.8.0]

### 🐛 **Bug Fix: Issue #90**

**Problem**: Nucleosome/MSP coordinates at alignment block boundaries were incorrectly extended into subsequent blocks, causing artificial length inflation.

### 🚀 **New Features: PanGenome Commands**

#### 1. `pg-inject`

Create mock BAM files from reference FASTA with perfectly aligned sequences while injecting annotation information from BED files.

#### 3. `pg-pansn`

Add/strip panSN-spec prefixes from BAM contig names

- Full panSN specification compliance
- Bidirectional prefix management

#### 2. `pg-lift`

Lift annotations through pangenome graphs from source to target coordinates

- Coordinate transformation across pangenome references
- Maintains annotation integrity during lifting
- Uses `vg` + `pg-pansn` + `pg-lift` under the hood.

### 📦 **New Core Architecture**

#### Major New Modules

- **`fibertig.rs`** (647 lines) - Core pangenome annotation framework
- **`bamannotations.rs`** (573 lines) - BED annotation handling & coordinate mapping
- **`panspec.rs`** (63 lines) - PanSN specification support

#### Enhanced Modules

- **`bamlift.rs`** - Improved coordinate lifting with better error handling
- **`bio_io.rs`** - Enhanced I/O capabilities (+114 lines)
- **`fiber.rs`** - Major cleanup and simplification (-218 lines)

#### Removed Legacy Code

- **`bamranges.rs`** (355 lines removed) - Replaced with more efficient implementations

#### Python (py-ft) Enhancements

- fiberplot - Standalone CLI tool for generating Altair HTML visualizations
- Updated dependencies and improved compatibility across the Python package

## [0.7.0] - 2025-07-10

### 🚀 Features

- Add `ft validate` subcommand for Fiber-seq BAM validation
- Add `-u` flag to all commands to allow uncompressed BAM output, this is very useful for piping Fiber-seq data between tools.
- Add haplotype support to footprint command
- Enable Revio model for Vega sequencing platform
- Allow uppercase letters for the mod code in modBAM parsing

### 🔧 Improvements

- Upgrade to burn 0.18 for improved ML model performance
- Outline Fiber-HMM functionality (experimental)
- Extensive clippy fixes and code cleanup

## [0.6.4] - 2025-01-31

- Fix: Error on unrecognized PacBio chemistry.

## [0.6.3] - 2025-01-31

- Fix: bit flag filtering #70.

## [0.6.1] - 2024-12-11

- Fix: the cli version string.

## [0.6.0] - 2024-11-14

- fix: do not report m6A predictions that happen within the first 7 or last 7 bp of a read. This is so the ML model only operates on real data. No changes to other calls. Will fix #65
- fix: report footprint codes even if there is no spanning msp, fixes #63
- feat: add a pyft utility that can take extract --all data and make it long format for plotting.
- feat: Add a shuffle option to pileup to help with the FDR calculations in FIRE
- feat: make nucleosomes and MSPs optional in pileup
- chore: use vergen for cli version
- feat: add phasing stats to QC
- feat: allow strip-base mods to filter on base mod quality.

## [0.5.4] - 2024-09-10

### ⚙️ Miscellaneous Tasks

- Update changelog with release
- Readme
- Readme
- Readme
- Release fibertools-rs version 0.5.4

## [0.5.3] - 2024-07-31

### ⚙️ Miscellaneous Tasks

- Update changelog with release
- Update changelog with release
- Update changelog with release
- Update changelog with release
- Update changelog with release
- Release fibertools-rs version 0.5.3

## [0.5.2] - 2024-07-31

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.5.2

## [0.5.1] - 2024-07-31

### 🐛 Bug Fixes

- Remove windows target
- Citation

### Release

- Remove windows target

## [0.4.2] - 2024-03-21

### ⚙️ Miscellaneous Tasks

- Clean
- Clean
- Clean cargo config
- Update helps
- Clean changelog
- Clean docs
- Clean docs

## [0.4.1] - 2024-03-20

### 🚀 Features

- Convert pytorch models to onnx to allow more backends (#46)

### ⚙️ Miscellaneous Tasks

- Sync lock with branch
- Sync lock with branch
- Ignore notebook
- Make test data smaller
- Make cargo publish smaller
- Make cargo publish smaller
- Make cargo publish smaller
- Make cargo publish smaller
- Make cargo publish smaller
- Make test data larger again since I remove the whole dir from cargo publish
- Make test data larger again since I remove the whole dir from cargo publish
- Make test data larger again since I remove the whole dir from cargo publish
- Make test data larger again since I remove the whole dir from cargo publish

## [0.3.9] - 2024-02-28

### 🚀 Features

- Add optional min msp filter to fire output
- Add optional min msp filter to fire output
- Add optional min msp filter to fire output
- Add optional min msp filter to fire output

### 🐛 Bug Fixes

- Allow overwritting fire qual
- Typo

### ⚙️ Miscellaneous Tasks

- Pyft lock
- Pyft lock
- Pyft lock
- Tch update
- Tch update
- Release fibertools-rs version 0.3.9

## [0.3.8] - 2024-01-24

### 🚀 Features

- Update pyft
- Update pyft to have a writer function
- Update pyft to have a writer function
- Update pyft to have a writer function
- Better fire paralization
- Center and extract now include fire scores in the range of 0-255

### 🐛 Bug Fixes

- GBDT made a semvar breaking change, fixing and locking depandancies.

### ⚙️ Miscellaneous Tasks

- Improved error msg
- Clippy
- CI
- CI
- Simplify colos
- Whitespace
- Readme and help pages
- Readme and help pages
- Pyft lock

## [0.3.7] - 2024-01-10

### 🚀 Features

- Add hp tag to center
- Add hp tag to center
- Add simplify options to center
- Move prediction into a feature only avalible through Cargo, since bioocnda wont let me do large installs anymore.
- Add ability to print FIRE features to a text file, TODO predict FIREs in rust.
- Fire progress indicator.
- Add a min ML options to add-nucs and keep working on FIRE predictions
- Speed up writing FIRE results
- FIRE io
- Start working on decorators
- Start working on decorators
- Use threads better, fire feats
- Use threads better, fire feats
- Use threads better, fire feats
- Use threads better, fire feats
- Update fire model
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- Adding rle information to fire feats
- FIRE now works in fibertools! TODO add info to extract, etc.
- Adding a fire bed+ extract for the fire pipeline.
- FIRE extract update
- FIRE extract update

### 🐛 Bug Fixes

- Try and better use threads.
- Try and better use threads.
- Make an iterator for fiberseq records
- Allow for the unaligned read start and read end
- Allow for the unaligned read start and read end

### ⚙️ Miscellaneous Tasks

- Update.
- Update.
- Update.
- Update.
- Update.
- Update.
- Update.
- Update.
- Add trace to ml mm parser

### Feat

- Adding code outline for footprinting tool

### Chroe

- Clippy
- Cli

## [0.3.6] - 2023-10-09

### ⚙️ Miscellaneous Tasks

- Update

## [bio-io-v0.3.2] - 2023-10-09

### ⚙️ Miscellaneous Tasks

- Improve progress bar
- Improve progress bar

## [bio-io-v0.3.1] - 2023-10-08

### Chrom

- Add docs.

## [0.3.3] - 2023-10-08

### 🐛 Bug Fixes

- Progress bar display issues.

### Chrom

- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.
- Add docs.

## [0.3.2] - 2023-09-27

### 🐛 Bug Fixes

- Extra comma in ft center wide format. fixes 21

### ⚙️ Miscellaneous Tasks

- Fix typo in cli message
- Update deps

### Chrom

- Fix warning msg

## [bio-io-v0.3.0] - 2023-08-13

### 🚀 Features

- Better progress bar for extract.
- Move progress bar into bamchunk iterator, unify progress bar.
- Move progress bar into bamchunk iterator, unify progress bar.
- Move progress bar into bamchunk iterator, unify progress bar.
- Move progress bar into bamchunk iterator, unify progress bar.

## [bio-io-v0.2.0] - 2023-07-29

### ⚙️ Miscellaneous Tasks

- Release bio-io version 0.2.0

## [bamlift-v0.2.0] - 2023-07-29

### 🚀 Features

- Very large improvment in speed of lifting over ranges, some liftover results differ from before by 1bp but it is rare.
- Adding liftover to pyft
- Reorganize pyft api and docs
- Adding ft center to the python module!
- Add qual to ft cetner, and clean the ft center code more
- Change the result of liftovers to be an Option<i64>.

### 🐛 Bug Fixes

- Move nuc and msp logic out of extract into ranges struct and simplify. We ran on a whole genome ft extract to confirm results dont change.
- Reduce number of copies for nuc and msp.
- Unify api anming for liftover, optimize for inclusion in pyft

### ⚙️ Miscellaneous Tasks

- Example py-ft
- Hide more deps under features for pyft
- Clippy
- Clippy
- Readme
- Readme
- Readme
- Readme
- Readme
- Readme
- Readme
- Readme
- Readme
- Readme
- Regex update
- Fix docs and fetch
- Drop darkmode
- Improve docs
- Refactor a large part of the code base to reduce redudance. TODO simplify center and basemods with this new api.
- Simplifed center with new api.
- Release bamlift version 0.2.0

## [0.2.6] - 2023-07-21

### 🚀 Features

- I have a minimal working python package yay!
- Add cpg to python modele and make easier access in the rust basemod api

### 🐛 Bug Fixes

- Rework fiberdata to consume bam record to avoid extra copy
- Speed options

### ⚙️ Miscellaneous Tasks

- Add to changelog
- Add to changelog
- Readme
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Python docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Docs
- Rename iterator
- Rename iterator
- Rename iterator
- Rename iterator

## [0.2.5] - 2023-07-18

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.5

## [0.2.4] - 2023-07-11

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.4

## [0.2.3] - 2023-06-24

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.3

## [0.2.2] - 2023-06-15

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.2

## [0.2.1] - 2023-06-14

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.1

## [0.2.0] - 2023-06-10

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.2.0

## [0.1.4] - 2023-04-19

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.4

## [0.1.3] - 2023-03-14

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.3

## [0.1.2] - 2023-02-13

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.2

## [0.1.1] - 2023-02-02

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.1-alpha.4
- Release fibertools-rs version 0.1.1

## [0.1.1-alpha.3] - 2023-02-01

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.1-alpha.3

## [0.1.1-alpha.2] - 2023-02-01

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.1-alpha.1
- Release fibertools-rs version 0.1.1-alpha.2

## [0.1.1-alpha] - 2023-02-01

### ⚙️ Miscellaneous Tasks

- Release fibertools-rs version 0.1.1-alpha

<!-- generated by git-cliff -->
