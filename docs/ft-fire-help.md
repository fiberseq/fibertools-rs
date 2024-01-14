```
Add FIREs (Fiber-seq Inferred Regulatory Elements) to a bam file with m6a predictions

Usage: ft fire [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Bam HiFi file with m6A and MSP calls [default: -]
  [OUT]  Output file (bam by default, table if --feats_to_text is used, and bed9 + if --extract is used) [default: -]

Options:
  -e, --extract                                                                        Output just FIRE elements in bed9 format
  -s, --skip-no-m6a                                                                    Don't write reads with no m6A calls to the output bam
  -w, --width-bin <WIDTH_BIN>                                                          Width of bin for feature collection [default: 40]
  -b, --bin-num <BIN_NUM>                                                              Number of bins to collect [default: 9]
      --best-window-size <BEST_WINDOW_SIZE>                                            Calculate stats for the highest X bp window within each MSP Should be a fair amount higher than the expected linker length [default: 100]
  -u, --use-5mc                                                                        Use 5mC data in FIREs
  -m, --min-msp-length-for-positive-fire-call <MIN_MSP_LENGTH_FOR_POSITIVE_FIRE_CALL>  Minium length of msp to call a FIRE [default: 85]
      --model <MODEL>                                                                  optional path to a model json file
      --fdr-table <FDR_TABLE>                                                          Optional path to a FDR table
  -f, --feats-to-text                                                                  Output FIREs features for training in a table format
  -h, --help                                                                           Print help
  -V, --version                                                                        Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```
