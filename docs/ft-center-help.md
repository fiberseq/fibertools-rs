```
This command centers fiberseq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs

Usage: ft center [OPTIONS] <BAM> <BED>

Arguments:
  <BAM>
          Fiberseq bam file, must be aligned and have an index

  <BED>
          Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file and corrected for strand if a 4th strand column is included

          If you include strand information in the 4th column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.

Options:
  -m, --min-ml-score <MIN_ML_SCORE>
          Minium score in the ML tag to include in the output

          [default: 150]

  -w, --wide
          Provide data in wide format, one row per read

  -h, --help
          Print help information (use `-h` for a summary)

  -V, --version
          Print version information

Global-Options:
  -t, --threads <THREADS>
          Threads for decompression

          [default: 8]

Debug-Options:
  -v, --verbose...
          Logging level [-v: Info, -vv: Debug, -vvv: Trace]

      --quiet
          Turn of all logging
```
