#!/usr/bin/env bash

for subcommand in "" "extract" "center"; do
    echo $subcommand
    out="docs/ft-${subcommand}-help.md"
    echo '```' >$out
    cargo run -- $subcommand --help >>$out
    echo '```' >>$out
done
