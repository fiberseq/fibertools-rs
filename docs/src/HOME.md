
<p align="center">
<img src="images/fiber_tools_teal.png" alt="fibertools-rs dark logo" width="200" class="center"/>
</p>

[![Actions Status](https://github.com/fiberseq/fibertools-rs/workflows/CI/badge.svg)](https://github.com/fiberseq/fibertools-rs/actions)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
[![Downloads](https://img.shields.io/conda/dn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)
[![crates.io version](https://img.shields.io/crates/v/fibertools-rs)](https://crates.io/crates/fibertools-rs)
[![crates.io downloads](https://img.shields.io/crates/d/fibertools-rs?color=orange&label=downloads)](https://crates.io/crates/fibertools-rs)
[![DOI](https://zenodo.org/badge/517338593.svg)](https://zenodo.org/badge/latestdoi/517338593)

This is the book for `fibertools` (`ft`) which is a CLI tool for **creating and interacting with Fiber-seq BAM** files. 

**Key features** include:
* [Predicting m6A](creatings/predict.md) sites from PacBio Fiber-seq data
* [Identifying FIREs](creating/fire.md) (**<ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements**)
* [Extracting](extracting/extract.md) Fiber-seq results into plain text files.
* [Centering](extracting/center.md) Fiber-seq results around a given position.
* [pyft](pyft.md): Python bindings for `fibertools`

A quick start guide can be found [here](quick-start.md), and a complete list of subcommands and their help pages can be found [here](help.md).

