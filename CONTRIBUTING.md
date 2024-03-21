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
```bash
cargo release --workspace {release type} -x
```
Where release type is one of: major, minor, patch.

I have also started trying cargo smart-release and it seems to work well:
```bash
cargo smart-release --update-crates-index --execute
```
