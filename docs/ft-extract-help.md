```
Extract fiberseq data into plain text files.

See https://fiberseq.github.io/fibertools-rs/docs/extract.html for a description of the outputs.

Usage: ft extract [OPTIONS] [BAM]

Arguments:
  [BAM]
          Input BAM file. If no path is provided extract will read bam data from stdin.
          
          For m6A prediction, this should be a HiFi bam file with kinetics data.
          
          For other commands, this should be a bam file with m6A calls.
          
          [default: -]

Options:
  -r, --reference
          Report in reference sequence coordinates

      --molecular
          Report positions in the molecular sequence coordinates

      --m6a <M6A>
          Output path for m6a bed12

  -c, --cpg <CPG>
          Output path for 5mC (CpG, primrose) bed12

      --msp <MSP>
          Output path for methylation sensitive patch (msp) bed12

  -n, --nuc <NUC>
          Output path for nucleosome bed12

  -a, --all <ALL>
          Output path for a tabular format including "all" fiberseq information in the bam

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

All-Format-Options:
  -q, --quality
          Include per base quality scores in "fiber_qual"

  -s, --simplify
          Simplify output by remove fiber sequence
```
