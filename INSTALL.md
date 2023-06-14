## First install `libtorch`
Get `libtorch` **v2.0.1** from the [PyTorch website](https://pytorch.org/get-started/) and extract the content of the zip file.
- On Linux/Unix system you can download with:
    * ```wget  https://download.pytorch.org/libtorch/cu116/libtorch-shared-with-deps-2.0.1%2Bcu118.zip```
- On macOS you can download with:
    * ```wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-2.0.1.zip```
- Windows is not supported and will not be.

Then add the following to your `.bashrc` or equivalent, where `/path/to/libtorch` is the path to the directory that was created when unzipping the file:
```bash
export LIBTORCH_CXX11_ABI=0
export LIBTORCH=/path/to/libtorch # e.g. export LIBTORCH=/Users/mrvollger/lib/libtorch
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```

## From `crates.io` [![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
Installation from `crates.io` requires the rust package manager `cargo`. You can find [how to install `cargo` here.](https://doc.rust-lang.org/cargo/getting-started/installation.html)
Furthermore, a recent version of `gcc` and `cmake` is required. I have tested and recommend `gcc v10.2.0` and `cmake v3.21.1`, though other versions may work.
```
cargo install fibertools-rs
```

## From `github` (active development)
```
cargo install --git https://github.com/mrvollger/fibertools-rs
```

