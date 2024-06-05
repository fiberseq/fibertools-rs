# Installation

This chapter covers the various ways to install `fibertools`.

It is easiest to install `fibertools` using conda ([instructions](install.md#from-bioconda)). However, for the fastest m6A predictions you need to install with `libtorch` enabled ([instructions](install.md#with-libtorch)) if you would like to install from source, you can use the following instructions.

# From `bioconda` [![Conda](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)

`fibertools-rs` is avalible through `bioconda` and can be installed with the following command:

```bash
mamba install -c conda-forge -c bioconda fibertools-rs
```

However, due to size constraints in `bioconda` this version does not support contain the pytorch libraries or GPU acceleration for m6A predictions. m6A predictions will still work in the bioconda version but may be much slower. If you would like to use m6A prediction and GPU acceleration, you will need to install using the directions in the [INSTALL.md](/INSTALL.md) file.


# From `crates.io` [![crates.io](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)

Installation from `crates.io` requires the rust package manager `cargo`. You can find [how to install `cargo` here.](https://doc.rust-lang.org/cargo/getting-started/installation.html). Furthermore, a recent version of `gcc` and `cmake` is required. I have tested and recommend `gcc v10.2.0` and `cmake v3.21.1`, though other versions may work.

```
cargo install fibertools-rs
```

# From GitHub (active development)
Using cargo from source:
```
cargo install --git https://github.com/fiberseq/fibertools-rs
```
or using git form source:
```bash
git clone https://github.com/fiberseq/fibertools-rs
cd fibertools-rs
cargo build --release
./target/release/ft --help
```



# With libtorch

Get `libtorch` **v2.2.0** from the [PyTorch website](https://pytorch.org/get-started/) and extract the content of the zip file.

- On Linux/Unix system you can download with:
  - `wget  https://download.pytorch.org/libtorch/cu118/libtorch-shared-with-deps-2.2.0%2Bcu118.zip`
- On macOS you can download with:
  - `wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-2.2.0.zip`
- Windows is not supported and will not be.

Then add the following to your `.bashrc` or equivalent, where `/path/to/libtorch` is the path to the directory that was created when unzipping the file:

```bash
export LIBTORCH_CXX11_ABI=0
export LIBTORCH=/path/to/libtorch # e.g. export LIBTORCH=/Users/mrvollger/lib/libtorch
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```

Finally install using `cargo`:
```bash
cargo install --all-features fibertools-rs
```