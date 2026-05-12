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

The changelog is managed by git-cliff which will run with the release action.

You can preview the changelog with:

```bash
git cliff | less
```

To cut a release, run:

```bash
cargo release {release type}
```

or

```bash
git commit -am "release: 0.2.0"
git tag "v0.2.0"
```

Where release type is one of: major, minor, patch.

The release page on GitHub will then be updated using `cargo dist`. You can preview this with:

```bash
cargo dist plan
```

and then run:

```bash
git push
git push --tags
```

### Re-triggering a failed release

If the release workflow fails, fix the issue (e.g. update `cargo-dist-version` in `dist-workspace.toml`, then run `dist generate-ci`), commit, push to main, and re-tag:

```bash
gh release delete v0.X.Y --yes
git push origin :refs/tags/v0.X.Y
git tag -d v0.X.Y
git tag v0.X.Y
git push origin v0.X.Y --force
```
