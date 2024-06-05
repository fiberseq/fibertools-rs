# From [![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/fibertools-rs?color=green)](https://anaconda.org/bioconda/fibertools-rs)

`fibertools-rs` is avalible through `bioconda` and can be installed with the following command:

```bash
mamba install -c conda-forge -c bioconda fibertools-rs
```

However, due to size constraints in `bioconda` this version does not support contain the pytorch libraries or GPU acceleration for m6A predictions. m6A predictions will still work in the bioconda version but may be much slower. If you would like to use m6A prediction and GPU acceleration, you will need to install using the directions in the [INSTALL.md](/INSTALL.md) file.