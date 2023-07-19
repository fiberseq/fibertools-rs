```
Strip out select base modifications

Usage: ft strip-basemods [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Bam HiFi file with base mods [default: -]
  [OUT]  Output bam file [default: -]

Options:
  -b, --basemod <BASEMOD>  base modification to strip out of the bam file [default: m6A] [possible
                           values: m6A, 6mA, 5mC, CpG]
  -h, --help               Print help
  -V, --version            Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
