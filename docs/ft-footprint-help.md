```
Infer footprints from fiberseq data

Usage: ft footprint [OPTIONS] <BAM> <BED> <YAML>

Arguments:
  <BAM>   Indexed and aligned bam file with m6A and MSP calls
  <BED>   BED file with the regions to footprint. Should all contain the same motif with proper strand information, and ideally be ChIP-seq peaks
  <YAML>  yaml describing the modules of the footprint

Options:
  -o, --out <OUT>  Output bam [default: -]
  -h, --help       Print help
  -V, --version    Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn off all logging
```
