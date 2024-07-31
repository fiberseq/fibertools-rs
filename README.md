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

`fibertools-rs` a CLI tool for creating and interacting with Fiber-seq BAM files. For more details read the [book](https://fiberseq.github.io/).

# Install [![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)

`fibertools-rs` is avalible through `bioconda` and can be installed with the following command:

```bash
mamba install -c conda-forge -c bioconda fibertools-rs
```

However, due to size constraints in `bioconda` this version does not support contain the pytorch libraries or GPU acceleration for m6A predictions. m6A predictions will still work in the bioconda version but may be much slower. If you would like to use m6A prediction and GPU acceleration, you will need to install using the directions in the [INSTALL.md](/INSTALL.md) file.

# Usage

```bash
ft --help
```

[Help page for fibertools](/docs/help.md)

# Highlighted subcommands for `fibertools-rs`

### `ft predict-m6a`

Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags. [Help page for predict-m6a](/docs/help.md#ft-predict-m6a).

### `ft add-nucleosomes`

Add nucleosomes to a bam that file already contains m6a predictions. Note, this process is also run in the background during `predict-m6a`, so it is unnecessary to run independently unless you want to try new parameters for nucleosome calling. [Help page for add-nucleosomes](/docs/help.md#ft-add-nucleosomes).

### `ft extract`

Extracts Fiber-seq data from a bam file into plain text. [Docs for extract](/docs/extract.md).

### `ft center`

Center Fiber-seq reads (bam) around reference position(s). [Docs for center](/docs/center.md).

### `ft footprint`
Footprint Fiber-seq reads (bam) around reference motifs(s). [Docs for footprint](/docs/footprint.md).

# Python API (`pyft`)

The python API is still in development and not stable; however, you can find the current code progress in the [py-ft](/py-ft) folder. More information available at [readthedocs](https://py-ft.readthedocs.io/en/latest/).

# Cite

**Jha, A.**, **Bohaczuk, S. C.**, Mao, Y., Ranchalis, J., Mallory, B. J., Min, A. T., Hamm, M. O., Swanson, E., Finkbeiner, C., Li, T., Whittington, D., Stergachis, A. B., & **Vollger, M. R.** (2023). DNA-m6A calling and integrated long-read epigenetic and genetic analysis with **fibertools**. _bioRxiv_. https://doi.org/10.1101/2023.04.20.537673

# Contributing
If you would like to contribute to `fibertools-rs`, please see the [CONTRIBUTING.md](/CONTRIBUTING.md) file for more information.
