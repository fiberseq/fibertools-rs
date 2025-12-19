# Development

## Setup

Create and activate a virtual environment in the python directory:

```bash
cd python
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Install maturin and dev dependencies:

```bash
pip install maturin pytest pysam
```

## Building

Build and install the extension in development mode:

```bash
maturin develop
```

After making changes to the Rust code, run `maturin develop` again to rebuild.

## Testing

Run the test suite:

```bash
pytest tests/ -v
```

## Cleaning

To clean build artifacts:

```bash
cargo clean
rm -f molecular_annotation/*.so
```

To fully reset (including uninstalling the package):

```bash
pip uninstall -y molecular-annotation
cargo clean
rm -f molecular_annotation/*.so
```
