```
Add FIREs (Fiber-seq Inferred Regulatory Elements) to a bam file with m6a predictions

Usage: ft fire [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]
          Input BAM file. If no path is provided extract will read bam data from stdin.
          
          For m6A prediction, this should be a HiFi bam file with kinetics data.
          
          For other commands, this should be a bam file with m6A calls.
          
          [default: -]

  [OUT]
          Output file (bam by default, table if --feats_to_text is used, and bed9 + if --extract is used)
          
          [default: -]

Options:
  -e, --extract
          Output just FIRE elements in bed9 format

  -s, --skip-no-m6a
          Don't write reads with no m6A calls to the output bam

      --min-msp <MIN_MSP>
          Skip reads without at least `N` MSP calls
          
          [env: MIN_MSP=]
          [default: 0]

      --min-ave-msp-size <MIN_AVE_MSP_SIZE>
          Skip reads without an average MSP size greater than `N`
          
          [env: MIN_AVE_MSP_SIZE=]
          [default: 0]

  -w, --width-bin <WIDTH_BIN>
          Width of bin for feature collection
          
          [env: WIDTH_BIN=]
          [default: 40]

  -b, --bin-num <BIN_NUM>
          Number of bins to collect
          
          [env: BIN_NUM=]
          [default: 9]

      --best-window-size <BEST_WINDOW_SIZE>
          Calculate stats for the highest X bp window within each MSP Should be a fair amount higher than the expected linker length
          
          [env: BEST_WINDOW_SIZE=]
          [default: 100]

      --min-msp-length-for-positive-fire-call <MIN_MSP_LENGTH_FOR_POSITIVE_FIRE_CALL>
          Minium length of msp to call a FIRE
          
          [env: MIN_MSP_LENGTH_FOR_POSITIVE_FIRE_CALL=]
          [default: 85]

      --model <MODEL>
          Optional path to a model json file. If not provided ft will use the default model (recommended)
          
          [env: FIRE_MODEL=]

      --fdr-table <FDR_TABLE>
          Optional path to a FDR table
          
          [env: FDR_TABLE=]

  -f, --feats-to-text
          Output FIREs features for training in a table format

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
