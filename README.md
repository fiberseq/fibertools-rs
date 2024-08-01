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
conda install bioconda::fibertools-rs
```

However, due to size constraints in `bioconda` this version does not support contain the pytorch libraries or GPU acceleration for m6A predictions. m6A predictions will still work in the bioconda version but may be much slower. If you would like to use m6A prediction and GPU acceleration, you will need to install using the directions in the fibertools [book](https://fiberseq.github.io/fibertools/install.html).

# Usage

```bash
ft --help
```

[Help page for fibertools](https://fiberseq.github.io/fibertools/help.html#ft)

# Highlighted subcommands for `fibertools-rs`

### `ft predict-m6a`

Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags. [Docs for predict-m6a](https://fiberseq.github.io/fibertools/creating/predict.html).

### `ft add-nucleosomes`

Add nucleosomes to a bam that file already contains m6a predictions. Note, this process is also run in the background during `predict-m6a`, so it is unnecessary to run independently unless you want to try new parameters for nucleosome calling. [Docs for add-nucleosomes](https://fiberseq.github.io/fibertools/creating/nucleosome.html).

### `ft extract`

Extracts Fiber-seq data from a bam file into plain text. [Docs for extract](https://fiberseq.github.io/fibertools/extracting/extract.html).

### `ft center`

Center Fiber-seq reads (bam) around reference position(s). [Docs for center](https://fiberseq.github.io/fibertools/extracting/center.html).

### `ft footprint`
Footprint Fiber-seq reads (bam) around reference motifs(s). [Docs for footprint](https://fiberseq.github.io/fibertools/extracting/footprint.html).

# Python API (`pyft`)

The python API is still in development and not stable; however, you can find the current code progress in the [py-ft](/py-ft) folder. More information available at [readthedocs](https://py-ft.readthedocs.io/en/latest/).

# Cite

**Jha, A.**, **Bohaczuk, S. C.**, Mao, Y., Ranchalis, J., Mallory, B. J., Min, A. T., Hamm, M. O., Swanson, E., Dubocanin, D., Finkbeiner, C., Li, T., Whittington, D., Noble, W. S., Stergachis, A. B., & **Vollger, M. R**. (2024). DNA-m6A calling and integrated long-read epigenetic and genetic analysis with fibertools. Genome Research. https://doi.org/10.1101/gr.279095.124

# Contributing
If you would like to contribute to `fibertools-rs`, please see the [CONTRIBUTING.md](/CONTRIBUTING.md) file for more information.
