```
Extract fiberseq data into plain text files

Usage: ft extract [OPTIONS] [BAM]

Arguments:
  [BAM]  Fiberseq bam file [default: -]

Options:
  -r, --reference                    Report in reference sequence coordinates
      --molecular                    Report positions in the molecular sequence coordinates
  -m, --min-ml-score <MIN_ML_SCORE>  Minium score in the ML tag to include in the output [default:
                                     125]
      --m6a <M6A>                    Output path for m6a bed12
  -c, --cpg <CPG>                    Output path for 5mC (CpG, primrose) bed12
      --msp <MSP>                    Output path for methylation sensitive patch (msp) bed12
  -n, --nuc <NUC>                    Output path for nucleosome bed12
  -a, --all <ALL>                    Output path for a tabular format including "all" fiberseq
                                     information in the bam
  -h, --help                         Print help
  -V, --version                      Print version

All-Format-Options:
  -q, --quality     Include per base quality scores in "fiber_qual"
  -f, --full-float  Add the full floating point predictions of the ML model
  -s, --simplify    Simplify output by remove fiber sequence

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
