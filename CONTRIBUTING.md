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

```bash
cargo release {release type} 
```
Where release type is one of: major, minor, patch.
