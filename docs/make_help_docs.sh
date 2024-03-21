#!/usr/bin/env bash

echo $LIBTORCH
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH

for subcommand in "" "extract" "center" "predict-m6a" "fire" "add-nucleosomes" "footprint" "clear-kinetics" "strip-basemods" "track-decorators"; do
    echo $subcommand
    out="docs/ft-${subcommand}-help.md"
    echo '```' >$out
    cargo run -- $subcommand --help >>$out
    echo '```' >>$out
done
