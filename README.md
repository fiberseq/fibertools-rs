# `fibertools-rs`
[![Actions Status](https://github.com/mrvollger/fibertools-rs/workflows/CI/badge.svg)](https://github.com/mrvollger/fibertools-rs/actions)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
 [![Downloads](https://img.shields.io/conda/dn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
[![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
[![crates.io downloads](https://img.shields.io/crates/d/fibertools-rs?color=orange&label=downloads)](https://crates.io/crates/fibertools-rs)
[![DOI](https://zenodo.org/badge/517338593.svg)](https://zenodo.org/badge/latestdoi/517338593)

`fibertools-rs` a CLI tool for creating and interacting with fiberseq bam files.

# Install
## From `crates.io` [![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
Installation from `crates.io` requires the rust package manager `cargo`. You can find [how to install `cargo` here.](https://doc.rust-lang.org/cargo/getting-started/installation.html)
Furthermore, a recent version of `gcc` and `cmake` is required. I have tested and recommend `gcc v10.2.0` and `cmake v3.21.1`, though other versions may work.
```
cargo install fibertools-rs
```

## From `bioconda` [![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
```
mamba install -c conda-forge -c bioconda fibertools-rs
```
## From `github` (active development)
```
cargo install --git https://github.com/mrvollger/fibertools-rs
```

# Usage
```bash
ft --help
```
[Help page for fibertools](/docs/ft--help.md)

# Subcommands for `fibertools-rs`
## `ft predict-m6a`
[Help page for predict-m6a](/docs/ft-predict-m6a-help.md). Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags.
### Installing with support for CNN m6A prediction
To allow for m6A predictions with the CNN model you must follow these modified installation instructions.
* Get `libtorch` **v1.12.0** from the [PyTorch website](https://pytorch.org/get-started/) download section and extract the content of the zip file.
    * On my linux system with a cuda gpu this is what I downloaded:
    * ```wget https://download.pytorch.org/libtorch/cu116/libtorch-cxx11-abi-shared-with-deps-1.12.0%2Bcu116.zip```
* Add the following to your `.bashrc` or equivalent, where `/path/to/libtorch` is the path to the directory that was created when unzipping the file:
```bash
export LIBTORCH=/path/to/libtorch # e.g. export LIBTORCH=/Users/mrvollger/lib/libtorch
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```
And install `fibertools-rs` from `cargo` with the `cnn` feature enabled:
```bash
cargo install fibertools-rs --features cnn
```
## `ft extract`
[Help page for extract](/docs/ft-extract-help.md). Extracts fiberseq data from a bam file into plain text.
### `ft extract --all`
The extract all option is a special option that tries to extract all the fiberseq data into a tabular format. The following is an image of the output. Note that the [column names](/docs/ft-all-columns.md) will be preserved across different software versions (unless otherwise noted); however, the order may change and new columns may be added. Therefore, when loading the data (with `pandas` e.g.) be sure to use the column names as opposed to indexes for manipulation.
![ft-extract all](/assets/img/ft-extract-all.png)


## `ft center`
[Help page for center](/docs/ft-center-help.md). Center fiberseq reads (bam) around reference position(s).
![Center](/assets/img/center.png)


# Read the fibertools docs
You can find the docs for the latest release here:
[https://docs.rs/fibertools-rs/latest/fibertools_rs/](https://docs.rs/fibertools-rs/latest/fibertools_rs/)
or download from source and run:
```
cargo doc --open --features cnn
```
and the docs will open in your browser.

# TODO for v0.0.11
- [ ] Add `rustybam` stats to ft `all` as an option
- [ ] add option result to bamlift
- [ ] Add more test cases, learn about test modules in folders
- [ ] Test GPU support, see if I can simplify or statically link PyTorch.
- [ ] Improve progress bar for predict-m6a.
    - [ ] Get size of bam, say how far we are through the bam in terms of MB/GB?
- [ ] Add unaligned, secondary, supplemental reads to the test bam.
- [ ] Detect GPU memory to set batch size dynamically.
