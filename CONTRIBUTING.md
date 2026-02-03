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
