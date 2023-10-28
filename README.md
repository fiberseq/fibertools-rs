---
---

# `fibertools-rs`

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

However, due to size constraints in `bioconda` this version does not support m6a prediction or GPU acceleration. If you would like to use m6A prediction and GPU acceleration, you will need to install using the directions in the [INSTALL.md](/INSTALL.md) file.

# Usage

```bash
ft --help
```

[Help page for fibertools](/docs/ft--help.md)

# Subcommands for `fibertools-rs`

## `ft predict-m6a`

Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags. [Help page for predict-m6a](/docs/ft-predict-m6a-help.md).

## `ft add-nucleosomes`

Add nucleosomes to a bam that file already contains m6a predictions. Note, this process is also run in the background during `predict-m6a`, so it is unnecessary to run independently unless you want to try new parameters for nucleosome calling. [Help page for add-nucleosomes](/docs/ft-add-nucleosomes-help.md).

## `ft extract`

Extracts Fiber-seq data from a bam file into plain text. [Help page for extract](/docs/extract.md).
![Extract](/assets/img/ft-extract-all.png)

## `ft center`

Center Fiber-seq reads (bam) around reference position(s). [Help page for center](/docs/center.md).
![Center](/assets/img/center.png)

# Python API (`pyft`)

The python API is still in development and not stable; however, you can find the current code progress in the [py-ft](/py-ft) folder. More information available at [readthedocs](https://py-ft.readthedocs.io/en/latest/).

# Cite

**Jha, A.**, **Bohaczuk, S. C.**, Mao, Y., Ranchalis, J., Mallory, B. J., Min, A. T., Hamm, M. O., Swanson, E., Finkbeiner, C., Li, T., Whittington, D., Stergachis, A. B., & **Vollger, M. R.** (2023). Fibertools: fast and accurate DNA-m6A calling using single-molecule long-read sequencing. _bioRxiv_. https://doi.org/10.1101/2023.04.20.537673

# Read the fibertools library docs

You can find the docs for the latest release here:
[https://docs.rs/fibertools-rs/latest/fibertools_rs/](https://docs.rs/fibertools-rs/latest/fibertools_rs/)
or download from source and run:

```
cargo doc --open
```

and the docs will open in your browser.

# TODO items

- [ ] Use new iterator for `ft extract` and group writes to try and improve the speed
- [ ] long format extract command
- [ ] Improve progress bar for predict-m6a.
  - [ ] Get size of bam, say how far we are through the bam in terms of MB/GB?
- [x] Add a python API (see py-ft for progress)
  - [x] extract api
  - [x] center api
  - [ ] improve docs
  - [ ] add default data viz
  - [ ] add conversion to pandas data frame or maybe anndata
- [x] GPU support
  - [ ] see if I can simplify or statically link PyTorch to get it onto bioconda
  - [ ] Detect GPU memory to set batch size dynamically.
- [ ] Add unaligned, secondary, supplemental reads to the test bam.
- [x] add option result to bamlift
- [x] improve speed of liftover closest in bamlift. It takes about 50% of the time.
- [x] Add more test cases, learn about test modules in folders
- [x] Set filters for ML depending on the model used

# Contributing

If you would like to contribute to `fibertools-rs`, please see the [CONTRIBUTING.md](/CONTRIBUTING.md) file for more information.
