```
Strip out select base modifications

Usage: ft strip-basemods [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]
          Input BAM file. If no path is provided extract will read bam data from stdin.
          
          For m6A prediction, this should be a HiFi bam file with kinetics data.
          
          For other commands, this should be a bam file with m6A calls.
          
          [default: -]

  [OUT]
          Output bam file
          
          [default: -]

Options:
  -b, --basemod <BASEMOD>
          base modification to strip out of the bam file
          
          [default: m6A]
          [possible values: m6A, 6mA, 5mC, CpG]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

BAM-Options:
  -F, --filter <BIT_FLAG>
          BAM bit flags to filter on, equivalent to `-F` in samtools view
          
          [default: 0]

      --ml <MIN_ML_SCORE>
          Minium score in the ML tag to use or include in the output
          
          [default: 125]

Global-Options:
  -t, --threads <THREADS>
          Threads
          
          [default: 8]

Debug-Options:
  -v, --verbose...
          Logging level [-v: Info, -vv: Debug, -vvv: Trace]

      --quiet
          Turn off all logging
```
