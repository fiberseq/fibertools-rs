```
Make decorated bed files for fiberseq data

Usage: ft track-decorators [OPTIONS] --bed12 <BED12> [BAM]

Arguments:
  [BAM]  Bam HiFi file with m6A calls [default: -]

Options:
  -b, --bed12 <BED12>                Output path for bed12 file to be decorated
  -d, --decorator <DECORATOR>        Output path for decorator bed file [default: -]
  -m, --min-ml-score <MIN_ML_SCORE>  Minium score in the ML tag to include in the output [default: 125]
  -h, --help                         Print help
  -V, --version                      Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```
