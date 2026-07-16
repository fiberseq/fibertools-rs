# Migrate release automation: release-please → release-plz

**Date:** 2026-07-15
**Status:** Plan for review (do not execute until approved)

## Why

release-please's `release-type: "rust"` **locksteps the Cargo workspace** — on
every release it upgrades *all* workspace members to the release version. It
bumped `molecular-annotation` to `0.10.1` (PR #113) while leaving the
`[workspace.dependencies]` pin at `^0.0.1`, breaking `cargo` resolution and
`main`'s CI. This is fundamentally incompatible with our decoupled design
(MA at `0.0.x`, fibertools at `0.10.x`).

[release-plz](https://release-plz.dev) is Rust/Cargo-native: it versions each
crate **independently** based on that crate's own changed commits, updates
internal dependency requirements correctly, is semver-aware (optional
`cargo-semver-checks`), and **publishes to crates.io itself**. It would not have
produced this incident.

## Already done (recovery step 1)

- `main` unbroken: `molecular-annotation/Cargo.toml` restored to `0.0.1` to match
  the pin (commit `8ea2880`). No `v0.10.1` tag or crates.io release was ever
  created, so nothing was published in error.

## Validated locally with the release-plz CLI (0.3.160)

Ran `release-plz update` on a scratch branch (baseline: `fibertools-rs 0.10.0`,
`molecular-annotation 0.0.1`, matching crates.io). Findings:

- ✅ **Decoupled versioning works.** The mock_fire `fix:` (fibertools-only)
  bumped **only** fibertools `0.10.0 → 0.10.1`; MA stayed `0.0.1`. No lockstep.
- ✅ **Pin updated correctly.** When MA does bump, release-plz rewrites
  `[workspace.dependencies] molecular-annotation = { …, version = "…" }` — the
  exact thing release-please broke.
- ⚠️ **Prerequisite 1 — backfill baseline tags.** release-plz uses
  `<pkg>-v<version>` git tags as each crate's baseline. `v0.10.0` exists, but
  `molecular-annotation-v0.0.1` was **missing** (the manual `cargo publish`
  never tagged it), so release-plz spuriously bumped MA to `0.0.2` — counting
  the release-please version-churn commits as MA changes. Creating
  `molecular-annotation-v0.0.1` fixed it: MA correctly stayed put.
- ⚠️ **Prerequisite 2 — .gitignore hygiene.** release-plz refuses to run while
  files are **both committed and gitignored**. Offenders:
  `molecular-annotation/test-data/**` (from `molecular-annotation/.gitignore`),
  `py-ft/Cargo.lock`, `Train-FIRE/pixi.toml`, `pixi.lock` (root `.gitignore`).
  Must untrack them or remove the ignore rules before the Action will work.

## The load-bearing decision (answer first)

**Do we keep cargo-dist prebuilt binaries + the `curl | sh` installer?** This is
the *only* thing that adds complexity — release-plz publishes to crates.io
inline with just `GITHUB_TOKEN`.

- **Option S — Drop cargo-dist (simplest).** One `release-plz.yml` workflow,
  `GITHUB_TOKEN` only, no PAT, no dispatch, no draft handshake. release-plz opens
  the release PR, and on merge bumps versions, updates pins, tags, publishes to
  crates.io, and cuts GitHub releases. Users install via `cargo install` or
  `cargo binstall`. **Loses prebuilt binaries + shell installer.**
- **Option B — Keep cargo-dist binaries.** Same release-plz core, plus we keep
  `dispatch-releases = true` and have release-plz's release step `gh workflow
  run` cargo-dist (no PAT), with release-plz creating a **draft** release that
  cargo-dist attaches to and un-drafts. More moving parts (the same handshake
  that's bitten us), needs live validation.

The plan below covers **Option B** (keep binaries) and notes where **Option S**
simply omits steps. Recommendation: if binaries aren't heavily used, take
Option S — it removes every failure mode we've hit.

## Secondary decisions

1. **Token.** release-plz core needs only `GITHUB_TOKEN`. Tradeoff with the
   default token: CI checks won't run *on the release PR*, and tag/release
   events won't trigger other workflows. Since release-plz publishes crates.io
   inline and (Option B) triggers cargo-dist via `workflow_dispatch` (which
   `GITHUB_TOKEN` may fire), **no PAT is required**. Accept "no CI on the release
   PR", or add a PAT/App later if wanted.
2. **crates.io auth.** Reuse Trusted Publishing: mint an OIDC token with
   `rust-lang/crates-io-auth-action` and pass it as `CARGO_REGISTRY_TOKEN` to the
   release-plz release job. (Both crates are already published and TP-configured
   for MA; add TP for fibertools-rs.)
3. **Version baseline.** `main` currently has `fibertools-rs = 0.10.1` (from the
   #113 bump) but crates.io has `0.10.0`. Reset `main`'s `fibertools-rs` to
   `0.10.0` so release-plz starts from the published baseline and cleanly
   proposes `0.10.1` from the mock_fire `fix:` commit — rather than guessing
   about a half-applied bump.

---

## Phase 0 — Tear down release-please + prerequisites

- Delete `.github/workflows/release-please.yml`, `release-please-config.json`,
  `.release-please-manifest.json`.
- **Backfill the MA baseline tag** (Prerequisite 1): create and push
  `molecular-annotation-v0.0.1` pointing at the commit whose MA source matches
  the published `0.0.1` crate. (`v0.10.0` already exists for fibertools.)
- **Fix .gitignore hygiene** (Prerequisite 2): for each committed-and-ignored
  file, either `git rm --cached` it (if it shouldn't be tracked) or remove the
  matching ignore rule (test fixtures like `molecular-annotation/test-data/**`
  are intentionally committed → un-ignore those). Confirm
  `git ls-files -ci --exclude-standard` is empty afterward.
- Delete the stray `release-please--branches--main` branch if it reappears, and
  confirm no `release-please` PR is open (the stuck "untagged merged release PR"
  state disappears once the workflow/config are gone).
- Reset `main` versions to the published baseline: `fibertools-rs = 0.10.0`
  (root `Cargo.toml`), `molecular-annotation = 0.0.1` (already), pin `^0.0.1`
  (already). Update `Cargo.lock`.
- Keep `dist-workspace.toml` as-is (Option B) — it already has
  `dispatch-releases = true`, `create-release = false`, and MA `dist = false`.
  (Option S: delete cargo-dist entirely — `dist-workspace.toml`, the generated
  `release.yml`.)

## Phase 1 — release-plz config (`release-plz.toml`)

Independent per-crate versioning is release-plz's default (it is NOT lockstep,
unlike release-please's rust strategy). Per-package overrides only:

```toml
[[package]]
name = "fibertools-rs"
git_tag_name = "v{{ version }}"          # keep v0.X.Y (matches existing tags + cargo-dist)
git_release_draft = true                 # draft: cargo-dist attaches binaries, then un-drafts
features_always_increment_minor = true   # feat -> minor even on 0.x

[[package]]
name = "molecular-annotation"
git_tag_name = "molecular-annotation-v{{ version }}"
# git_release_draft defaults false -> MA release is published immediately
#   (it ships no binaries, so nothing un-drafts it).
# features_always_increment_minor defaults false -> keeps MA low (feat -> patch).

# semver_check is left OFF for now: it needs cargo-semver-checks installed in CI,
# and conventional-commit bumping is sufficient. Enable later by adding a
# cargo-semver-checks install step and `semver_check = true`.
```

- The `git_tag_name` values line up with the Phase 0 baselines: `v0.10.0`
  (exists) and `molecular-annotation-v0.0.1` (backfilled), so release-plz finds
  each crate's last release correctly.
- Validated by the CLI dry-runs: fibertools-only change → only fibertools bumps;
  MA change → MA bumps, fibertools re-releases, pin updated.

## Phase 2 — release-plz workflow (`.github/workflows/release-plz.yml`)

One workflow, `GITHUB_TOKEN` only (no PAT). release-plz publishes to crates.io
inline; the single cross-workflow hop (to cargo-dist) uses `workflow_dispatch`,
which `GITHUB_TOKEN` may trigger.

```yaml
name: release-plz

on:
  push:
    branches: [main]

permissions:
  contents: write
  pull-requests: write

jobs:
  # Open/update the release PR (version bumps + changelog + pin updates).
  release-plz-pr:
    runs-on: ubuntu-latest
    concurrency:
      group: release-plz-${{ github.ref }}
      cancel-in-progress: false
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0            # release-plz needs full history
          submodules: recursive
      - uses: dtolnay/rust-toolchain@stable
      - uses: release-plz/action@v0.5
        with:
          command: release-pr
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}

  # On merge of the release PR: publish unpublished crates to crates.io
  # (dependency-ordered), push tags, cut releases; then fire cargo-dist for
  # fibertools-rs.
  release-plz-release:
    runs-on: ubuntu-latest
    permissions:
      contents: write
      id-token: write     # crates.io Trusted Publishing (OIDC)
      actions: write      # dispatch cargo-dist (release.yml)
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0
          submodules: recursive
      - uses: dtolnay/rust-toolchain@stable
      - name: Authenticate to crates.io (Trusted Publishing)
        id: auth
        uses: rust-lang/crates-io-auth-action@v1
      - name: release-plz release
        id: release-plz
        uses: release-plz/action@v0.5
        with:
          command: release
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          CARGO_REGISTRY_TOKEN: ${{ steps.auth.outputs.token }}
      # Build binaries + installer for fibertools-rs, attach to its DRAFT
      # release, and un-draft it. Only runs when fibertools-rs was released.
      - name: Dispatch cargo-dist for fibertools-rs
        if: ${{ steps.release-plz.outputs.releases_created == 'true' }}
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          RELEASES: ${{ steps.release-plz.outputs.releases }}
        run: |
          tag=$(echo "$RELEASES" | jq -r '.[] | select(.package_name=="fibertools-rs") | .tag')
          if [ -n "$tag" ] && [ "$tag" != "null" ]; then
            echo "Dispatching cargo-dist for $tag"
            gh workflow run release.yml --ref main -f tag="$tag"
          else
            echo "fibertools-rs not in this release; skipping cargo-dist"
          fi
```

Notes / tradeoffs:
- **`GITHUB_TOKEN` only:** CI won't run *on the release PR*, and release-plz's
  tag push won't trigger tag-triggered workflows. Neither matters here —
  publishing is inline and cargo-dist is dispatched explicitly. Add a PAT/App
  later only if you want CI on the release PR.
- **crates.io Trusted Publishing:** `crates-io-auth-action` mints the OIDC token
  → `CARGO_REGISTRY_TOKEN`. Needs TP configured on crates.io for **both** crates
  (MA done; add fibertools-rs).
- **Draft handshake:** release-plz pushes the real git tag itself (not GitHub's
  lazy tag), so cargo-dist can check it out — no `force-tag-creation` needed
  (unlike release-please). `git_release_draft = true` (fibertools) gives
  cargo-dist a draft to attach to and un-draft.

## Phase 3 — Validation (live, on a scratch/pre-release)

1. release-plz opens a release PR proposing **only** `fibertools-rs 0.10.1`
   (MA untouched), with the pin correct — confirms no lockstep.
2. Merge it: crates.io gets `fibertools-rs 0.10.1`; tag + GitHub release created.
3. Option B: cargo-dist binaries attach to the release and it un-drafts.
4. Make a throwaway MA-only change → release-plz bumps only MA (e.g. `0.0.2`)
   and updates fibertools's pin — confirms decoupled versioning end-to-end.

## What this fixes vs the old setup

- No workspace lockstep → MA stays `0.0.x`, decoupled, as intended.
- Dependency pins updated automatically and correctly.
- Semver-aware bumps → no accidental `1.0.0`.
- crates.io publish is native → deletes the gated-job publish workflow.
- Still no PAT required.

## Open questions for you

1. **Option S or B** — drop cargo-dist binaries (simplest), or keep them?
2. **Token** — `GITHUB_TOKEN` only (no CI on the release PR), acceptable?
3. **Baseline reset** — OK to set `main`'s `fibertools-rs` back to `0.10.0` so
   release-plz cleanly proposes `0.10.1`?
