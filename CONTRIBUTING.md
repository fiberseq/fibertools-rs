# Contributing to `fibertools-rs`

Please feel free to open PRs! But first make sure you code passes tests, and please add tests for new features:

```bash
cargo test --all-features
```

Also format your code and check it with clippy before submitting a PR:

```bash
cargo fmt
cargo clippy --workspace
```

## Regression tests

End-to-end regression tests live in `tests/regression/` and compare subcommand
output against committed snapshots using [insta](https://insta.rs/). Each test
runs the `ft` binary against a fixture in `tests/data/`, projects a stable
subset of columns (so trailing/added columns aren't false positives), and
asserts the result matches the snapshot in `tests/regression/snapshots/`.

Run just the regression suite with:

```bash
cargo test --test regression
```

When you intentionally change output, regenerate snapshots with
[`cargo-insta`](https://insta.rs/docs/cli/):

```bash
cargo install cargo-insta      # one-time
cargo insta test --review      # walk through diffs interactively
# or, if the new output is correct as-is:
cargo insta accept
```

Inspect the resulting `.snap` diff before committing.

## Cutting a release

Releases are automated with [release-plz](https://release-plz.dev), which is
Cargo-native and versions each crate **independently**. You do not bump versions
or create tags by hand.

The two crates version separately:

- `fibertools-rs` (the `ft` CLI) — `v0.X.Y` tags, ships prebuilt binaries.
- `molecular-annotation` (the library) — `molecular-annotation-v0.0.z` tags,
  kept at a low internal version.

`fibertools-rs` pins an exact `molecular-annotation` version, so a change to the
library bumps the library **and** re-releases `fibertools-rs` against it (with
the pin updated automatically). A change that only touches `fibertools-rs`
bumps only `fibertools-rs`.

### How a release happens

1. Merge PRs to `main` using [Conventional
   Commits](https://www.conventionalcommits.org/). On these `0.x` crates:
   `fix:` → patch, `feat:` → minor (fibertools) / patch (molecular-annotation),
   `feat!:`/`BREAKING CHANGE:` → minor (capped below `1.0.0`). Squash-merge is
   recommended so each PR is one logical change.
2. release-plz keeps a **release PR** open (`chore: release`) that bumps the
   changed crates, updates the changelogs, and updates the internal dependency
   pin. Which crate bumps is decided by which crate's files changed.
3. When you are ready to release, **merge the release PR**. release-plz then:
   - publishes the changed crates to crates.io (dependency-ordered) via Trusted
     Publishing;
   - pushes the git tags and cuts GitHub releases — the `fibertools-rs` release
     is created as a **draft**;
   - dispatches `cargo dist` (the `Release` workflow) for `fibertools-rs`, which
     builds the `ft` binaries + shell installer, uploads them to the draft
     release, and **un-drafts** it.

Changelogs (`CHANGELOG.md`, `molecular-annotation/CHANGELOG.md`) are generated
by release-plz — do not edit them by hand.

### Previewing a release (dry run) with the release-plz CLI

Before merging, you can see exactly what release-plz would do, locally:

```bash
pixi global install release-plz          # or: cargo binstall release-plz

# `release-plz update` rewrites versions + changelogs in your working tree.
# Inspect the proposal, then revert it — this is a preview, not a commit.
release-plz update
git --no-pager diff                      # what would change (versions, pins, changelogs)
git checkout -- Cargo.toml Cargo.lock molecular-annotation/   # revert
rm -f CHANGELOG.md molecular-annotation/CHANGELOG.md          # remove generated changelogs
```

For the preview to be accurate, three things must hold (they normally do on a
clean `main`):

- **Full git history** is present (a shallow clone will misreport versions).
- **Baseline tags exist** for the currently published versions: `v0.X.Y`
  (fibertools-rs) and `molecular-annotation-v0.0.z`. release-plz uses these to
  know each crate's last release; a missing tag makes it propose a spurious bump.
- **No file is both committed and gitignored** — release-plz refuses to run
  otherwise. Check with `git ls-files -ci --exclude-standard` (must be empty).

### Trusted Publishing

crates.io publishing uses Trusted Publishing (OIDC, no stored token). Both
crates must have a trusted publisher configured on crates.io for repo
`fiberseq/fibertools-rs`, workflow **`release-plz.yml`**.

### Re-triggering a failed release

- If **crates.io publish** partially failed, just re-run the `release-plz` job —
  `release-plz release` is idempotent and skips already-published versions.
- If the **`cargo dist` (`Release`) run** failed (or the `fibertools-rs` release
  is stuck as a draft), fix the cause and re-dispatch it from the Actions tab
  (Run workflow → `tag = v0.X.Y`). It uploads to the existing draft and
  un-drafts when the artifacts are attached.
