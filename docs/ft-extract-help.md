```
ft-extract 0.0.8
Extract fiberseq data into plain text files

USAGE:
    ft extract [OPTIONS] [BAM]

ARGS:
    <BAM>    fiberseq bam file [default: -]

OPTIONS:
    -r, --reference
            report in reference sequence coordinates

    -m, --min-ml-score <MIN_ML_SCORE>
            Minium score in the ML tag to include in the output [default: 150]

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
            Print help information

    -V, --version
            Print version information

ALL-FORMATTING-OPTIONS:
    -q, --quality       include per base quality scores in "fiber_qual"
    -f, --full-float    Add the full floating point predictions of the ML model
    -s, --simplify      simplify output by remove fiber sequence

GLOBAL-OPTIONS:
    -t, --threads <THREADS>    Threads for decompression [default: 8]

DEBUG-OPTIONS:
    -v, --verbose    Logging level [-v: Info, -vv: Debug, -vvv: Trace]
        --quiet      Turn of all logging
```
