# Examples

## fiberseq_to_ma

Converts fiber-seq CRAM tags (ns/nl, as/al/aq) to the Molecular Annotation format.

```bash
# Run with default test file
cargo run --example fiberseq_to_ma

# Run with custom input
cargo run --example fiberseq_to_ma -- path/to/file.cram
```

### Input tags (fiber-seq format)
- `ns` / `nl` - Nucleosome starts and lengths
- `as` / `al` / `aq` - MSP/accessible region starts, lengths, and qualities

### Output tags
- **MA** - Inline format with start-length pairs
- **AQ** - Quality scores
- **M2** - Alternative separate format (starts only)
- **AL** - Lengths array for M2

Both formats are written to enable compression comparison.