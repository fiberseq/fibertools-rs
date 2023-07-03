```
Fiber-seq toolkit in rust

Usage: ft [OPTIONS] <COMMAND>

Commands:
  predict-m6a      Predict m6A positions using HiFi kinetics data and encode the results in the
                       MM and ML bam tags. Also adds nucleosome (nl, ns) and MTase sensitive patches
                       (al, as) [aliases: m6A, m6a]
  add-nucleosomes  Add nucleosomes to a bam file with m6a predictions
  extract          Extract fiberseq data into plain text files [aliases: ex, e]
  center           This command centers fiberseq data around given reference positions. This is
                       useful for making aggregate m6A and CpG observations, as well as
                       visualization of SVs [aliases: c, ct]
  clear-kinetics   Remove HiFi kinetics tags from the input bam file
  strip-basemods   Strip out select base modifications
  help             Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Global-Options:
  -t, --threads <THREADS>  Threads [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
