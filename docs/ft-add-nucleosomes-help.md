```
Add nucleosomes to a bam file with m6a predictions

Usage: ft add-nucleosomes [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Bam HiFi file with m6A calls [default: -]
  [OUT]  Output bam file with nucleosome calls [default: -]

Options:
  -n, --nucleosome-length <NUCLEOSOME_LENGTH>
          Minium nucleosome length [default: 75]
  -c, --combined-nucleosome-length <COMBINED_NUCLEOSOME_LENGTH>
          Minium nucleosome length when combining over a single m6A [default: 100]
  -d, --distance-from-end <DISTANCE_FROM_END>
          Minimum distance from the end of a fiber to call a nucleosome or MSP [default: 45]
  -h, --help
          Print help information
  -V, --version
          Print version information

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
