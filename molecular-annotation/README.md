# molecular-annotation

> Specification: see [here](https://github.com/fiberseq/Molecular-annotation-spec) for the MolecularAnnotation tag format.

A Rust library for reading and writing Molecular Annotation (MA) tags in SAM/BAM/CRAM files.

## Installation

```toml
[dependencies]
molecular-annotation = { path = "." }

# For BAM/CRAM file support
molecular-annotation = { path = ".", features = ["htslib"] }
```

## Encoding Formats

Two encoding formats are available via feature flags:

- **`inline-lengths`** (default): Lengths in MA string: `MA:Z:1000;nuc+:100-50,200-60`
- **`separate-lengths`**: Lengths in AL array: `MA:Z:1000;nuc+:100,200` + `AL:B:I,50,60`

## Documentation

```bash
cargo doc --open
```

## Examples

See [examples/README.md](examples/README.md) for available examples.

## License

MIT
