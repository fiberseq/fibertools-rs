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
However, the `bioconda` version currently does not support GPU acceleration. If you would like to use GPU acceleration, you will need to install using the directions in the [INSTALL.md](/INSTALL.md) file.


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
Extracts fiberseq data from a bam file into plain text. [Help page for extract](/docs/ft-extract-help.md). 
### `ft extract --all`
The extract all option is a special option that tries to extract all the fiberseq data into a tabular format. The following is an image of the output. Note that the [column names](/docs/ft-all-columns.md) will be preserved across different software versions (unless otherwise noted); however, the order may change and new columns may be added. Therefore, when loading the data (with `pandas` e.g.) be sure to use the column names as opposed to indexes for manipulation.
![ft-extract all](/assets/img/ft-extract-all.png)
## `ft center`
Center fiberseq reads (bam) around reference position(s). [Help page for center](/docs/ft-center-help.md).
![Center](/assets/img/center.png)

# Cite
**Jha, A.**, **Bohaczuk, S. C.**, Mao, Y., Ranchalis, J., Mallory, B. J., Min, A. T., Hamm, M. O., Swanson, E., Finkbeiner, C., Li, T., Whittington, D., Stergachis, A. B., & **Vollger, M. R.** (2023). Fibertools: fast and accurate DNA-m6A calling using single-molecule long-read sequencing. *bioRxiv*. https://doi.org/10.1101/2023.04.20.537673

# Read the fibertools library docs
You can find the docs for the latest release here:
[https://docs.rs/fibertools-rs/latest/fibertools_rs/](https://docs.rs/fibertools-rs/latest/fibertools_rs/)
or download from source and run:
```
cargo doc --open
```
and the docs will open in your browser.

# TODO items
- [ ] Add a python API (see py-ft for progress)
- [ ] Use new iterator for `ft extract` and group writes to try and improve the speed
- [x] Set filters for ML depending on the model used
- [ ] long format extract command
- [ ] add option result to bamlift
- [ ] improve speed of liftover closest in bamlift. It takes about 50% of the time. 
- [x] Add more test cases, learn about test modules in folders
- [x] GPU support
    - [ ] see if I can simplify or statically link PyTorch to get it onto bioconda
    - [ ] Detect GPU memory to set batch size dynamically.
- [ ] Improve progress bar for predict-m6a.
    - [ ] Get size of bam, say how far we are through the bam in terms of MB/GB?
- [ ] Add unaligned, secondary, supplemental reads to the test bam.