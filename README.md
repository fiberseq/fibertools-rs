# `fibertools-rs`
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
 [![Downloads](https://img.shields.io/conda/dn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
[![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
[![crates.io downloads](https://img.shields.io/crates/d/fibertools-rs?color=orange&label=downloads)](https://crates.io/crates/fibertools-rs)
[![DOI](https://zenodo.org/badge/517338593.svg)](https://zenodo.org/badge/latestdoi/517338593)

`fibertools-rs` a CLI tool for interacting with fiberseq bam files.

# Install
## From [![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs) (recommended)
```
cargo install fibertools-rs
```
[How to install `cargo`.](https://doc.rust-lang.org/cargo/getting-started/installation.html)

## From [![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs) (latest version not always available)
```
mamba install -c bioconda fibertools-rs
```
## From `github` (active development)
```
cargo install --git https://github.com/mrvollger/fibertools-rs
```



# `ft extract`
Extracts fiberseq data from a bam file into plain text.
```bash
ft-extract 0.0.4
Extract fiberseq data into plain text files

USAGE:
    ft extract [OPTIONS] [BAM]

ARGS:
    <BAM>    fiberseq bam file [default: -]

OPTIONS:
    -r, --reference    report in reference sequence coordinates
        --m6a <M6A>    Output path for m6a bed12
    -c, --cpg <CPG>    Output path for 5mC (CpG, primrose) bed12
        --msp <MSP>    Output path for methylation sensitive patch (msp) bed12
    -n, --nuc <NUC>    Output path for nucleosome bed12
    -a, --all <ALL>    Output path for
    -h, --help         Print help information
    -V, --version      Print version information
```


# ft-center
Center a fiberseq reads (bam) around a reference position(s).
```bash
ft-center 0.0.4
This command centers fiberseq data around given reference positions. This is useful for making
aggregate m6a and CpG observations, as well as visualization of SVs

USAGE:
    ft center [OPTIONS] <BAM> <BED>

ARGS:
    <BAM>    fiberseq bam file, must be aligned and have an index
    <BED>    Bed file on which to center fiberseq reads. Data is adjusted to the start position
             of the bed file and corrected for strand if a 4th strand column is included

OPTIONS:
    -w, --wide       Provide data in wide format, one row per read
    -h, --help       Print help information
    -V, --version    Print version information
```
![center](/images/center.png)
