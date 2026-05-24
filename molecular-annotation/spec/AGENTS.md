# Code Agent Guide

This document helps AI agents get oriented when working on the Molecular Annotation Spec repository.

## Project Overview

This repository defines and implements the **MA (Molecular Annotation) tag specification** for SAM/BAM/CRAM files. It provides a standardized format for annotating non-genetic elements on individual DNA molecules (e.g., nucleosome positions, methylation-accessible patches, regulatory elements from Fiber-seq data).

### Repository Structure

```
Molecular-annotation-spec/
├── README.md              # The specification document
├── rust/                  # Core Rust library
│   ├── src/
│   │   ├── lib.rs         # Main MolecularAnnotations type and API
│   │   ├── types.rs       # Core types (Annotation, Strand, QualityType, etc.)
│   │   ├── liftover.rs    # Coordinate liftover between query and reference
│   │   └── tests.rs       # Unit tests
│   └── Cargo.toml
├── python/                # Python bindings (PyO3)
│   ├── src/lib.rs         # PyO3 wrapper around Rust library
│   ├── molecular_annotation/
│   │   ├── __init__.py    # Package exports
│   │   └── pysam_utils.py # Utilities for pysam integration
│   ├── tests/
│   └── pyproject.toml
└── test_data/             # Sample BAM files for testing
```

## Key Concepts

### Coordinate Systems

**Critical**: The library uses two coordinate systems:

| Context       | System                               | Example                                |
| ------------- | ------------------------------------ | -------------------------------------- |
| Internal API  | **0-based half-open** `[start, end)` | Position 100, length 50 → `[100, 150)` |
| MA tag string | **1-based closed**                   | Same annotation → `101-50` in tag      |

The `from_tags()`/`to_tags()` methods automatically convert between these.

### Coordinate Orientations

For reverse-aligned reads:

- `get_coords()` / `get_bam_coords()` - BAM orientation (flipped for reverse reads)
- `get_forward_coords()` / `get_molecular_coords()` - Original molecular orientation (never flipped)

### MA Tag Format

```
MA:Z:read_length;type1+P:start1-len1,start2-len2;type2-:start3-len3
```

- Read length comes first
- Annotation types: `name` + strand (`+`/`-`/`.`) + optional quality type (`P`/`Q`)
- Positions are 1-based, lengths follow after `-`

### Related Tags

- **AQ:B:C** - Quality scores (u8 array, only for types with `P` or `Q`)
- **AN:Z** - Annotation names (comma-separated, empty strings for unnamed)

## Development Commands

### Rust

```bash
cd rust
cargo test              # Run all tests
cargo test --features htslib  # With BAM integration
cargo clippy            # Lint
cargo doc --open        # Generate docs
```

### Python

**Important**: Follow the dev notes in `python/DEV.md`:

```bash
cd python
source .venv/bin/activate
env -u CONDA_PREFIX maturin develop  # Build and install
python -m pytest tests/ -v           # Run tests
```

## Architecture Notes

### Rust Library (`rust/src/`)

- **lib.rs**: `MolecularAnnotations` - main container type

  - Holds `read_length`, `annotation_types`, optional `aligned_blocks`
  - Methods for adding/getting annotations, serialization, liftover

- **types.rs**: Core data types

  - `Annotation` - single annotation (start, length, quality, name)
  - `AnnotationType` - group of annotations with same type/strand/quality_type
  - `Strand`, `QualityType`, `Encoding` enums
  - `ParseError` for error handling

- **liftover.rs**: Coordinate transformation
  - `AlignedBlocks` - holds alignment blocks from CIGAR
  - `lift_to_reference()` / `lift_to_query()` - bidirectional liftover
  - Handles gaps (insertions/deletions) with snapping behavior

### Python Bindings (`python/src/lib.rs`)

- PyO3 wrapper exposing the same API
- Adds Pythonic conveniences like `ends` parameter (alternative to `lengths`)
- `pysam_utils.py` - integration with pysam for reading/writing BAM records

## Common Tasks

### Adding a New Feature

1. Implement in Rust first (`rust/src/`)
2. Add tests to `rust/src/tests.rs`
3. Expose via Python bindings (`python/src/lib.rs`)
4. Add Python tests (`python/tests/`)
5. Update README.md spec if needed

### Modifying the Spec

1. Update `README.md` with new format/behavior
2. Update Rust parser/serializer in `lib.rs`
3. Update Python bindings
4. Ensure all tests pass

### Working with Coordinates

- All API methods use **0-based half-open** intervals
- The only place 1-based appears is in MA tag strings
- Use `u32` for all coordinates and lengths (changed from i64)
- Reference coordinates can be `Option<u32>` (None if outside aligned region)

## Testing

Current test counts:

- **77 Rust tests** (including liftover edge cases)
- **11 Rust doc-tests**
- **25 Python tests**

Run both before committing:

```bash
cd rust && cargo test
cd ../python && source .venv/bin/activate && python -m pytest tests/ -v
```

## Recent Changes (as of Feb 2025)

1. **Removed legacy AL tag** - Lengths are now inline in MA tag (`start-length` format)
2. **Converted i64 to u32** - All coordinates/lengths now use u32
3. **Added method aliases** - `get_bam_coords()`, `get_molecular_coords()`
4. **Improved coordinate documentation** - Clear 0-based vs 1-based distinction

## Tips for AI Agents

1. **Always read the file before editing** - The codebase has specific patterns
2. **Use the dev notes** - `python/DEV.md` has the correct build commands
3. **Run tests after changes** - Both Rust and Python
4. **Coordinate systems matter** - Double-check 0-based vs 1-based in context
5. **The spec is in README.md** - Keep implementation aligned with spec
