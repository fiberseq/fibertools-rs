```
Make a pileup track of Fiber-seq features from a FIRE bam

Usage: ft pileup [OPTIONS] [BAM] [RGN]

Arguments:
  [BAM]
          Input BAM file. If no path is provided extract will read bam data from stdin.
          
          For m6A prediction, this should be a HiFi bam file with kinetics data.
          
          For other commands, this should be a bam file with m6A calls.
          
          [default: -]

  [RGN]
          Region string to make a pileup of. If not provided will make a pileup of the whole genome

Options:
  -o, --out <OUT>
          Output file
          
          [default: -]

  -m, --m6a
          include m6A calls

  -c, --cpg
          include 5mC calls

      --haps
          For each column add two new columns with the hap1 and hap2 specific data

  -k, --keep-zeros
          Keep zero coverage regions

  -p, --per-base
          Write output one base at a time even if the values do not change

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
