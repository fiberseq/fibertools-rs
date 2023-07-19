#!/usr/bin/env bash

echo $LIBTORCH
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH

for subcommand in "" "extract" "center" "predict-m6a" "add-nucleosomes" "clear-kinetics" "strip-basemods"; do
    echo $subcommand
    out="docs/ft-${subcommand}-help.md"
    echo '```' >$out
    cargo run -- $subcommand --help >>$out
    echo '```' >>$out
done

# extract the column names for the all command
printf '# Columns in `ft extract --all`\n```\n' >docs/ft-all-columns.md
cargo run -- -v extract --all - tests/data/all.bam | head -n 1 | sed 's/\t/\n/g' >>docs/ft-all-columns.md
echo '```' >>docs/ft-all-columns.md
