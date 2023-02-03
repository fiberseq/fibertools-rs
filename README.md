---
---
`fibertools-rs`
==============

<img src="./assets/img/fiber_tools_teal.png#gh-dark-mode-only" alt="fibertools-rs dark logo" width="200"/>
<img src="./assets/img/fiber_tools_grey.png#gh-light-mode-only" alt="fibertools-rs light logo" width="200"/>


[![Actions Status](https://github.com/fiberseq/fibertools-rs/workflows/CI/badge.svg)](https://github.com/fiberseq/fibertools-rs/actions)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
 [![Downloads](https://img.shields.io/conda/dn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
[![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
[![crates.io downloads](https://img.shields.io/crates/d/fibertools-rs?color=orange&label=downloads)](https://crates.io/crates/fibertools-rs)
[![DOI](https://zenodo.org/badge/517338593.svg)](https://zenodo.org/badge/latestdoi/517338593)

`fibertools-rs` a CLI tool for creating and interacting with fiberseq bam files.

# Install [![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)

`fibertools-rs` is avalible through `bioconda` and can be installed with the following command:
```bash
mamba install -c conda-forge -c bioconda fibertools-rs
```
Other installation methods are available in the [INSTALL.md](/INSTALL.md) file.


# Usage
```bash
ft --help
```
[Help page for fibertools](/docs/ft--help.md)

# Subcommands for `fibertools-rs`
## `ft predict-m6a`
[Help page for predict-m6a](/docs/ft-predict-m6a-help.md). Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags.
### Adding nucleosome calls to the BAM files
To add nucleosome calls to the BAM files you can use the python package [fibertools](https://github.com/fiberseq/fibertools#add-nucleosomes-and-msps-to-a-fibertools-rs-m6a-bam). See that repository for installation and instructions.

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
cargo doc --open
```
and the docs will open in your browser.

# TODO items
- [ ] Use new iterator for `ft extract` and group writes to try and improve the speed
- [x] Set filters for ML depending on the model used
- [ ] long format extract command
- [ ] Add `rustybam` stats to ft `all` as an option
- [ ] add option result to bamlift
- [ ] Add more test cases, learn about test modules in folders
- [x] Test GPU support, see if I can simplify or statically link PyTorch.
- [ ] Improve progress bar for predict-m6a.
    - [ ] Get size of bam, say how far we are through the bam in terms of MB/GB?
- [ ] Add unaligned, secondary, supplemental reads to the test bam.
- [ ] Detect GPU memory to set batch size dynamically.
