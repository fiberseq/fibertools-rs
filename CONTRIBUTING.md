# Contributing to `fibertools-rs`
Please feel free to open PRs! But first make sure you code passes tests, and please add tests for new features:
```bash
cargo test -p bio-io -p bamlift -p fibertools-rs
```
Also format your code and check it with clingy:
```bash
cargo fmt 
cargo clippy --workspace
```

## Cutting a release
```bash
cargo release --workspace {release type} -x
```
Where release type is one of: major, minor, patch.