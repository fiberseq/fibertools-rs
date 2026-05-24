# Development

## Setup

Create and activate a virtual environment in the python directory:

```bash
cd python
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Ensure pip is available and install dev dependencies:

```bash
python -m ensurepip --upgrade
python -m pip install maturin pytest pysam
```

**Note:** Use `python -m pip` instead of just `pip` to ensure you're using the venv's pip, not a system/conda pip that may be earlier in PATH.

## Building

Build and install the extension in development mode:

```bash
maturin develop
```

After making changes to the Rust code, run `maturin develop` again to rebuild.

### Conda Environment Conflicts

If you have a conda environment active alongside the venv, you may encounter build errors due to `CONDA_PREFIX` interfering with maturin's build process. To resolve this, either:

1. Deactivate conda first:
   ```bash
   conda deactivate
   source .venv/bin/activate
   maturin develop
   ```

2. Or unset `CONDA_PREFIX` for the build command:
   ```bash
   env -u CONDA_PREFIX maturin develop
   ```

## Testing

Run the test suite:

```bash
python -m pytest tests/ -v
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
