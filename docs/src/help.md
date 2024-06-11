# Help pages for fibertools subcommands
<!-- toc -->

## `ft `
```console
Fiber-seq toolkit in rust

Usage: ft [OPTIONS] <COMMAND>

Commands:
  predict-m6a       Predict m6A positions using HiFi kinetics data and encode the results in the
                        MM and ML bam tags. Also adds nucleosome (nl, ns) and MTase sensitive
                        patches (al, as) [aliases: m6A, m6a]
  add-nucleosomes   Add nucleosomes to a bam file with m6a predictions
  fire              Add FIREs (Fiber-seq Inferred Regulatory Elements) to a bam file with m6a
                        predictions
  extract           Extract fiberseq data into plain text files [aliases: ex, e]
  center            This command centers fiberseq data around given reference positions. This is
                        useful for making aggregate m6A and CpG observations, as well as
                        visualization of SVs [aliases: c, ct]
  footprint         Infer footprints from fiberseq data
  track-decorators  Make decorated bed files for fiberseq data
  pileup            Make a pileup track of Fiber-seq features from a FIRE bam
  clear-kinetics    Remove HiFi kinetics tags from the input bam file
  strip-basemods    Strip out select base modifications
  help              Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft predict-m6a`
```console
Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags.
Also adds nucleosome (nl, ns) and MTase sensitive patches (al, as)

Usage: ft predict-m6a [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [OUT]  Output bam file with m6A calls in new/extended MM and ML bam tags [default: -]

Options:
  -n, --nucleosome-length <NUCLEOSOME_LENGTH>
          Minium nucleosome length [default: 75]
  -c, --combined-nucleosome-length <COMBINED_NUCLEOSOME_LENGTH>
          Minium nucleosome length when combining over a single m6A [default: 100]
      --min-distance-added <MIN_DISTANCE_ADDED>
          Minium distance needed to add to an already existing nuc by crossing an m6a [default: 25]
  -d, --distance-from-end <DISTANCE_FROM_END>
          Minimum distance from the end of a fiber to call a nucleosome or MSP [default: 45]
  -k, --keep
          Keep hifi kinetics data
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging

Developer-Options:
      --force-min-ml-score <FORCE_MIN_ML_SCORE>
          Force a different minimum ML score
      --all-calls
          Keep all m6A calls regardless of how low the ML value is
  -b, --batch-size <BATCH_SIZE>
          Number of reads to include in batch prediction [default: 1]
```

## `ft fire`
```console
Add FIREs (Fiber-seq Inferred Regulatory Elements) to a bam file with m6a predictions

Usage: ft fire [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [OUT]  Output file (BAM by default, table of MSP features if `--feats-to-text` is used, and bed9 +
         if `--extract`` is used) [default: -]

Options:
  -e, --extract
          Output just FIRE elements in bed9 format
      --all
          When extracting bed9 format include all MSPs and nucleosomes
  -f, --feats-to-text
          Output FIREs features for training in a table format
  -s, --skip-no-m6a
          Don't write reads with no m6A calls to the output bam
      --min-msp <MIN_MSP>
          Skip reads without at least `N` MSP calls [env: MIN_MSP=] [default: 0]
      --min-ave-msp-size <MIN_AVE_MSP_SIZE>
          Skip reads without an average MSP size greater than `N` [env: MIN_AVE_MSP_SIZE=] [default:
          0]
  -w, --width-bin <WIDTH_BIN>
          Width of bin for feature collection [env: WIDTH_BIN=] [default: 40]
  -b, --bin-num <BIN_NUM>
          Number of bins to collect [env: BIN_NUM=] [default: 9]
      --best-window-size <BEST_WINDOW_SIZE>
          Calculate stats for the highest X bp window within each MSP Should be a fair amount higher
          than the expected linker length [env: BEST_WINDOW_SIZE=] [default: 100]
      --min-msp-length-for-positive-fire-call <MIN_MSP_LENGTH_FOR_POSITIVE_FIRE_CALL>
          Minium length of msp to call a FIRE [env: MIN_MSP_LENGTH_FOR_POSITIVE_FIRE_CALL=]
          [default: 85]
      --model <MODEL>
          Optional path to a model json file. If not provided ft will use the default model
          (recommended) [env: FIRE_MODEL=]
      --fdr-table <FDR_TABLE>
          Optional path to a FDR table [env: FDR_TABLE=]
  -h, --help
          Print help
  -V, --version
          Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft extract`
```console
Extract fiberseq data into plain text files

Usage: ft extract [OPTIONS] [BAM]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]

Options:
  -r, --reference  Report in reference sequence coordinates
      --molecular  Report positions in the molecular sequence coordinates
      --m6a <M6A>  Output path for m6a bed12
  -c, --cpg <CPG>  Output path for 5mC (CpG, primrose) bed12
      --msp <MSP>  Output path for methylation sensitive patch (msp) bed12
  -n, --nuc <NUC>  Output path for nucleosome bed12
  -a, --all <ALL>  Output path for a tabular format including "all" fiberseq information in the bam
  -h, --help       Print help (see more with '--help')
  -V, --version    Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging

All-Format-Options:
  -q, --quality   Include per base quality scores in "fiber_qual"
  -s, --simplify  Simplify output by remove fiber sequence
```

## `ft center`
```console
This command centers fiberseq data around given reference positions. This is useful for making
aggregate m6A and CpG observations, as well as visualization of SVs

Usage: ft center [OPTIONS] --bed <BED> [BAM]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]

Options:
  -b, --bed <BED>    Bed file on which to center fiberseq reads. Data is adjusted to the start
                     position of the bed file and corrected for strand if the strand is indicated in
                     the 6th column of the bed file. The 4th column will also be checked for the
                     strand but only after the 6th is. If you include strand information in the 4th
                     (or 6th) column it will orient data accordingly and use the end position of bed
                     record instead of the start if on the minus strand. This means that profiles of
                     motifs in both the forward and minus orientation will align to the same central
                     position
  -d, --dist <DIST>  Set a maximum distance from the start of the motif to keep a feature
  -w, --wide         Provide data in wide format, one row per read
  -r, --reference    Return relative reference position instead of relative molecular position
  -s, --simplify     Replace the sequence output column with just "N"
  -h, --help         Print help (see more with '--help')
  -V, --version      Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft add-nucleosomes`
```console
Add nucleosomes to a bam file with m6a predictions

Usage: ft add-nucleosomes [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [OUT]  Output bam file with nucleosome calls [default: -]

Options:
  -n, --nucleosome-length <NUCLEOSOME_LENGTH>
          Minium nucleosome length [default: 75]
  -c, --combined-nucleosome-length <COMBINED_NUCLEOSOME_LENGTH>
          Minium nucleosome length when combining over a single m6A [default: 100]
      --min-distance-added <MIN_DISTANCE_ADDED>
          Minium distance needed to add to an already existing nuc by crossing an m6a [default: 25]
  -d, --distance-from-end <DISTANCE_FROM_END>
          Minimum distance from the end of a fiber to call a nucleosome or MSP [default: 45]
  -h, --help
          Print help
  -V, --version
          Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft footprint`
```console
Infer footprints from fiberseq data

Usage: ft footprint [OPTIONS] --bed <BED> --yaml <YAML> [BAM]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]

Options:
  -b, --bed <BED>    BED file with the regions to footprint. Should all contain the same motif with
                     proper strand information, and ideally be ChIP-seq peaks
  -y, --yaml <YAML>  yaml describing the modules of the footprint
  -o, --out <OUT>    Output bam [default: -]
  -h, --help         Print help
  -V, --version      Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft clear-kinetics`
```console
Remove HiFi kinetics tags from the input bam file

Usage: ft clear-kinetics [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [OUT]  Output bam file without hifi kinetics [default: -]

Options:
  -h, --help     Print help
  -V, --version  Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft strip-basemods`
```console
Strip out select base modifications

Usage: ft strip-basemods [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [OUT]  Output bam file [default: -]

Options:
  -b, --basemod <BASEMOD>  base modification to strip out of the bam file [default: m6A] [possible
                           values: m6A, 6mA, 5mC, CpG]
  -h, --help               Print help
  -V, --version            Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft track-decorators`
```console
Make decorated bed files for fiberseq data

Usage: ft track-decorators [OPTIONS] --bed12 <BED12> [BAM]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]

Options:
  -b, --bed12 <BED12>          Output path for bed12 file to be decorated
  -d, --decorator <DECORATOR>  Output path for decorator bed file [default: -]
  -h, --help                   Print help
  -V, --version                Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

## `ft pileup`
```console
Make a pileup track of Fiber-seq features from a FIRE bam

Usage: ft pileup [OPTIONS] [BAM] [RGN]

Arguments:
  [BAM]  Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a
         HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A
         calls [default: -]
  [RGN]  Region string to make a pileup of. If not provided will make a pileup of the whole genome

Options:
  -o, --out <OUT>   Output file [default: -]
  -m, --m6a         include m6A calls
  -c, --cpg         include 5mC calls
      --haps        For each column add two new columns with the hap1 and hap2 specific data
  -k, --keep-zeros  Keep zero coverage regions
  -p, --per-base    Write output one base at a time even if the values do not change
  -h, --help        Print help
  -V, --version     Print version

BAM-Options:
  -F, --filter <BIT_FLAG>  BAM bit flags to filter on, equivalent to `-F` in samtools view [default:
                           0]
      --ml <MIN_ML_SCORE>  Minium score in the ML tag to use or include in the output [default: 125]

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```

