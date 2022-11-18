#!/usr/bin/env bash

echo $LIBTORCH
#echo $LD_LIBRARY_PATH
#echo $DYLD_LIBRARY_PATH
#export LIBTORCH=/Users/mrvollger/lib/libtorch_1.12.0
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH

for subcommand in "" "extract" "center" "predict-m6a"; do
    echo $subcommand
    out="docs/ft-${subcommand}-help.md"
    echo '```' >$out
    cargo run -- $subcommand --help >>$out
    echo '```' >>$out
done

# extract the column names for the all command
printf '# Columns in `ft extract --all`\n```\n' >docs/ft-all-columns.md
cargo run -- -v extract --all - .test/all.bam | head -n 1 | sed 's/\t/\n/g' >>docs/ft-all-columns.md
echo '```' >>docs/ft-all-columns.md
