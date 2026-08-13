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

Lengths are encoded inline in the Ma string (`Ma:Z:1000;nuc+:100-50,200-60`).
The retired separate-length encoding (`Al` array) is stripped on write and
never emitted. Readers accept the all-uppercase `MA`/`AQ`/`AN` spellings that
fibertools 0.10-0.12 emitted; writers emit only `Ma`/`Aq`/`An`
(samtools/hts-specs#862).

## Documentation

```bash
cargo doc --open
```

## Examples

See [examples/README.md](examples/README.md) for available examples.

## License

MIT
