# Burn Framework Upgrade: 0.12 → 0.17.1 Benchmarking

This directory contains benchmarking results and binaries for upgrading fibertools-rs from Burn 0.12 to 0.17.1, maintaining compatibility with rayon multithreading for base modification prediction.

## Summary

Successfully upgraded the burn deep learning framework from version 0.12 to 0.17.1 while preserving full rayon parallelization. The key technical challenge was that burn 0.17's models are no longer `Sync`, breaking compatibility with rayon's `par_iter_mut()`. This was solved by creating fresh `PredictOptions` instances per thread.

## Technical Solution

### Problem

- Burn 0.17.1 models contain `OnceCell<Tensor>` which doesn't implement `Sync`
- Rayon's `par_iter_mut()` requires `Sync` traits for shared data

## Key Changes Made

### Dependencies (Cargo.toml)

- `burn = { version = "0.17.1", features = ["candle"] }`
- `burn-import = { version = "0.17.1" }`
- `tch = { version = "0.19.0", optional = true }`

### Code Updates

- **src/m6a_burn/mod.rs**: Updated TensorData API for burn 0.17 compatibility
- **src/subcommands/predict_m6a.rs**: Implemented factory pattern for thread safety

## Benchmarking Results

### Performance Comparison

| Version     | Backend     | Mean Time       | Relative Performance |
| ----------- | ----------- | --------------- | -------------------- |
| Burn 0.12   | PyTorch 2.2 | 2.201s ± 0.029s | 1.05x                |
| Burn 0.17.1 | PyTorch 2.6 | 2.087s ± 0.040s | **1.00x (baseline)** |
| Burn 0.12   | Candle      | 1.335s ± 0.025s | 1.56x faster         |
| Burn 0.17.1 | Candle      | 1.258s ± 0.028s | **1.66x faster**     |

### Key Findings

- **PyTorch versions**: Nearly identical performance (1.05x difference within margin of error)
- **Candle backend**: 1.66x faster than PyTorch, improved from 1.56x in 0.12
- **Identical output**: All 187 test reads produce identical base modification predictions (MM/ML tags)

## Output Verification

Verified that both versions produce identical results:

### Wrapper Scripts

- `ft_0_12_pytorch_wrapper.sh` - Sets up pytorch2.2 environment
- `ft_0_17_1_pytorch_factory_wrapper.sh` - Sets up pytorch2.6 environment

## Environment Setup

### PyTorch 2.6 (Burn 0.17.1)

```bash
conda activate pytorch2.6
export LIBTORCH_USE_PYTORCH=1
unset LIBTORCH
cargo build --release --all-features
```

### PyTorch 2.2 (Burn 0.12)

```bash
conda activate pytorch2.2
export LIBTORCH_USE_PYTORCH=1
unset LIBTORCH
cargo build --release --all-features
```
