```
Fiber-seq toolkit in rust

Usage: ft [OPTIONS] <COMMAND>

Commands:
  extract      Extract fiberseq data into plain text files [aliases: ex, e]
  center       This command centers fiberseq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs [aliases: c, ct]
  predict-m6a  Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags [aliases: m6A, m6a]
  help         Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help information
  -V, --version  Print version information

Global-Options:
  -t, --threads <THREADS>  Threads for decompression [default: 8]

Debug-Options:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
