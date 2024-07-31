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

First preview changes to the changelog:
```bash
git cliff --bump | less -S
```
Then if you agree with them bump the version:
```
git cliff --bump | less -S
```

```bash
cargo release {release type} 
```
Where release type is one of: major, minor, patch.
