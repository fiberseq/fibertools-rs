# fibertools-rs

`fibertools-rs` a CLI tool for interacting with fiberseq bam files. The main utility is `ft extract` which extracts fiberseq data from a bam file into plain text.

```bash
$ ft extract --help
ft-extract 0.0.1
Extract fiberseq data into plain text files

USAGE:
    ft extract [OPTIONS] [BAM]

ARGS:
    <BAM>    fiberseq bam file [default: -]

OPTIONS:
    -r, --reference    report in reference sequence coordinates
    -m, --m6a <M6A>    Output path for m6a bed12
    -c, --cpg <CPG>    Output path for CpG (primrose) bed12
    -m, --msp <MSP>    Output path for methylation sensitive patch (msp) bed12
    -n, --nuc <NUC>    Output path for nucleosome bed12
    -h, --help         Print help information
    -V, --version      Print version information
```
