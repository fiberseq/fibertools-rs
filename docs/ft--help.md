```
fibertools-rs 0.0.4
Mitchell R. Vollger <mrvollger@gmail.com>
fiberseq toolkit in rust

USAGE:
    ft [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -t, --threads <THREADS>    Threads for decompression [default: 8]
    -q, --quiet                Turn of all logging
    -h, --help                 Print help information
    -V, --version              Print version information

DEBUG:
    -v, --verbose    Logging level [-v: Debug, -vv: Trace]

SUBCOMMANDS:
    extract    Extract fiberseq data into plain text files [aliases: ex, e]
    center     This command centers fiberseq data around given reference positions. This is
                   useful for making aggregate m6a and CpG observations, as well as visualization of
                   SVs [aliases: c, ct]
    help       Print this message or the help of the given subcommand(s)
```
